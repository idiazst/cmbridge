.fit_sieve_md <- function(y, d, x, z, control) {
  ctrl <- .merge_control(list(
    target_basis = "bs", instrument_basis = "bs",
    target_degree = 3L, instrument_degree = 3L,
    target_df = 6L, instrument_df = 9L,
    lambda = 1e-8, weight_ridge = 1e-8,
    link = "identity", solver_tolerance = 1e-10, solver_max_iter = 10000L
  ), control)

  ctrl$link <- match.arg(ctrl$link, c("identity", "log", "inverse_logit"))
  if (!is.finite(ctrl$lambda) || ctrl$lambda < 0) stop("lambda must be nonnegative.", call. = FALSE)
  if (!is.finite(ctrl$solver_tolerance) || ctrl$solver_tolerance <= 0 ||
      !is.finite(ctrl$solver_max_iter) || ctrl$solver_max_iter < 1)
    stop("Invalid linked-sieve solver controls.", call. = FALSE)

  active <- abs(d) > 0
  if (!any(active)) stop("diagonal is zero for every observation.", call. = FALSE)
  if (anyNA(x[active, , drop = FALSE])) stop("target may be missing only where diagonal is zero.", call. = FALSE)
  if (anyNA(z)) stop("instrument cannot be missing for sieve minimum distance.", call. = FALSE)

  if (ctrl$link != "identity" && (ctrl$target_basis != "cell" || ctrl$instrument_basis != "cell"))
    stop("Positive sieve links require target_basis and instrument_basis = cell.", call. = FALSE)
  xspec <- .fit_basis_spec(x[active, , drop = FALSE], ctrl$target_basis,
                           ctrl$target_degree, ctrl$target_df)
  zspec <- .fit_basis_spec(z, ctrl$instrument_basis,
                           ctrl$instrument_degree, ctrl$instrument_df)
  if (xspec$type == "cell" && zspec$type == "cell") {
    return(.fit_sieve_cells(y, d, x, z, active, xspec, zspec, ctrl))
  }
  H <- matrix(0, nrow(x), ncol(.eval_basis_spec(x[active, , drop = FALSE], xspec)))
  H[active, ] <- .eval_basis_spec(x[active, , drop = FALSE], xspec)
  Q <- .eval_basis_spec(z, zspec)

  n <- length(y)
  A <- crossprod(Q, d * H) / n
  cvec <- crossprod(Q, y) / n
  W <- .safe_inverse(crossprod(Q) / n + ctrl$weight_ridge * diag(ncol(Q)), ridge = ctrl$weight_ridge)
  P <- diag(ncol(H)); P[1L, 1L] <- 0
  lhs <- crossprod(A, W %*% A) + ctrl$lambda * P
  rhs <- crossprod(A, W %*% cvec)
  coef <- as.numeric(.safe_solve(lhs, rhs))

  pred_train <- as.numeric(H %*% coef)
  pred_train[!active] <- NA_real_
  residual <- y - d * ifelse(is.na(pred_train), 0, pred_train)
  moment <- as.numeric(crossprod(Q, residual) / n)

  list(
    coefficients = coef,
    target_spec = xspec,
    instrument_spec = zspec,
    target_design_cols = ncol(H),
    fitted = pred_train,
    residual = residual,
    moment_loss = as.numeric(crossprod(moment, W %*% moment)),
    tuning = ctrl,
    predict_fun = .prediction_closure(function(newx) {
      newx <- .as_matrix(newx)
      as.numeric(.eval_basis_spec(newx, xspec) %*% coef)
    }, list(xspec = xspec, coef = coef)),
    moment_weight = W,
    instrument_features = .prediction_closure(function(newz) .eval_basis_spec(newz, zspec),
                                              list(zspec = zspec))
  )
}

# Compute exactly the same normal equations as the dense joint-category
# design. Only the small cross-tabulation is dense, regardless of sample size.
.fit_sieve_cells <- function(y, d, x, z, active, xspec, zspec, ctrl) {
  n <- length(y)
  xi <- match(.cell_keys(x[active, , drop = FALSE]), xspec$levels)
  zi <- match(.cell_keys(z), zspec$levels)
  nx <- length(xspec$levels); nz <- length(zspec$levels)
  J <- as.matrix(Matrix::sparseMatrix(i = zi[active], j = xi, x = d[active],
                                    dims = c(nz, nx)))
  A <- rbind(c(sum(J), colSums(J)), cbind(rowSums(J), J)) / n
  group_sum <- function(value) {
    grouped <- rowsum(value, zi, reorder = TRUE)
    out <- numeric(nz)
    out[as.integer(rownames(grouped))] <- grouped[, 1L]
    out
  }
  cvec <- c(sum(y), group_sum(y)) / n
  counts <- tabulate(zi, nbins = nz) / n
  covariance <- rbind(c(1, counts), cbind(counts, diag(counts, nrow = nz)))
  W <- .safe_inverse(covariance + ctrl$weight_ridge * diag(nz + 1L),
                     ridge = ctrl$weight_ridge)
  P <- diag(nx + 1L); P[1L, 1L] <- 0
  solver <- NULL
  if (ctrl$link == "identity") {
    coef <- as.numeric(.safe_solve(crossprod(A, W %*% A) + ctrl$lambda * P,
                                   crossprod(A, W %*% cvec)))
  } else {
    # Optimize unrestricted real link coefficients. The range restriction
    # follows only from the inverse-expit (or exponential) parameterization.
    # Retain the original ridge penalty on FUNCTION values for comparability.
    Acell <- A[, -1L, drop = FALSE]
    center <- diag(nx) - matrix(1 / nx, nx, nx)
    S <- crossprod(Acell, W %*% Acell) + ctrl$lambda * center
    b <- as.numeric(crossprod(Acell, W %*% cvec))
    solver <- .sieve_link_solve(S, b, as.numeric(crossprod(cvec, W %*% cvec)),
                              ctrl$link, ctrl$solver_tolerance, ctrl$solver_max_iter)
    values <- solver$values
    common <- mean(values)
    coef <- c(common, values - common)
  }
  pred_train <- rep(NA_real_, n)
  pred_train[active] <- if (ctrl$link == "identity") coef[1L] + coef[1L + xi] else solver$values[xi]
  residual <- y - d * ifelse(is.na(pred_train), 0, pred_train)
  moment <- c(sum(residual), group_sum(residual)) / n
  link_coef <- link_intercept <- NULL
  if (ctrl$link != "identity") {
    values <- solver$values
    link_coef <- solver$parameters
    log_gap <- if (ctrl$link == "log") link_coef else -link_coef
    log_mean_gap <- max(log_gap) + log(mean(exp(log_gap - max(log_gap))))
    link_intercept <- if (ctrl$link == "log") log_mean_gap else -log_mean_gap
    predict_fun <- .prediction_closure(function(newx) {
      newx <- .as_matrix(newx)
      if (ncol(newx) != xspec$p) stop("new data have the wrong number of columns.", call. = FALSE)
      index <- match(.cell_keys(newx), xspec$levels)
      eta <- rep(link_intercept, nrow(newx))
      seen <- !is.na(index)
      eta[seen] <- link_coef[index[seen]]
      if (link == "log") exp(eta) else 1 + exp(-eta)
    }, list(xspec = xspec, link = ctrl$link, link_coef = link_coef, link_intercept = link_intercept))
  } else {
    predict_fun <- .prediction_closure(function(newx) {
      newx <- .as_matrix(newx)
      if (ncol(newx) != xspec$p) stop("new data have the wrong number of columns.", call. = FALSE)
      index <- match(.cell_keys(newx), xspec$levels)
      coef[1L] + ifelse(is.na(index), 0, coef[1L + index])
    }, list(xspec = xspec, coef = coef))
  }
  list(
    coefficients = coef, link_coefficients = link_coef, link_intercept = link_intercept,
    solver = solver, target_spec = xspec, instrument_spec = zspec,
    target_design_cols = nx + 1L, fitted = pred_train, residual = residual,
    moment_loss = as.numeric(crossprod(moment, W %*% moment)), tuning = ctrl,
    predict_fun = predict_fun,
    moment_weight = W,
    instrument_features = .prediction_closure(function(newz) {
      newz <- .as_matrix(newz)
      if (ncol(newz) != zspec$p) stop("new data have the wrong number of columns.", call. = FALSE)
      index <- match(.cell_keys(newz), zspec$levels)
      seen <- which(!is.na(index)); nr <- nrow(newz)
      Matrix::sparseMatrix(i = c(seq_len(nr), seen),
                           j = c(rep(1L, nr), 1L + index[seen]), x = 1,
                           dims = c(nr, nz + 1L))
    }, list(zspec = zspec, nz = nz))
  )
}

# Inverse-expit is evaluated as 1 + exp(-eta) for numerical stability.
# Every optimizer parameter is an unrestricted real number; no box bounds,
# active-set solve, transformed coefficient clipping, or prediction clipping.
.sieve_link_objective <- function(S, b, constant, link) {
  sign <- if (link == "inverse_logit") -1 else 1
  offset <- if (link == "inverse_logit") 1 else 0
  scale <- max(diag(S))
  if (!is.finite(scale) || scale <= 0) stop("Linked sieve has no identified cell direction.", call. = FALSE)
  value <- function(eta) offset + exp(sign * eta)
  origin <- as.numeric(.safe_solve(S, b))
  origin_gradient <- as.numeric(S %*% origin - b)
  reference_loss <- as.numeric(constant + crossprod(origin, S %*% origin) - 2 * crossprod(b, origin))
  fn <- function(eta) {
    v <- value(eta)
    if (any(!is.finite(v))) return(Inf)
    difference <- v - origin
    as.numeric(crossprod(difference, S %*% difference) +
      2 * crossprod(origin_gradient, difference)) / scale
  }
  gr <- function(eta) {
    gap <- exp(sign * eta)
    2 * sign * gap * as.numeric(S %*% (offset + gap) - b) / scale
  }
  list(fn = fn, gr = gr, value = value, scale = scale, offset = offset, sign = sign,
       reference_loss = reference_loss)
}

.sieve_link_solve <- function(S, b, constant, link, tolerance, max_iter) {
  S <- (S + t(S)) / 2
  problem <- .sieve_link_objective(S, b, constant, link)
  # These values only initialize the unrestricted coefficients; the optimizer
  # may move in either direction over the entire real parameter space.
  raw <- as.numeric(.safe_solve(S, b))
  start <- problem$sign * log(pmax(.1, raw - problem$offset))
  initial_loss <- problem$fn(start)
  fits <- list()
  trace <- data.frame(method = character(), convergence = integer(), loss = double(),
                      coefficient_gradient = double(), function_gradient = double())
  assess <- function(par) {
    values <- problem$value(par)
    gradient <- 2 * as.numeric(S %*% values - b) / sqrt(diag(S))
    gap <- exp(problem$sign * par)
    # Near the limiting link value, inspect the function gradient as well as
    # the coefficient gradient: a nearly flat link must not hide a descent.
    function_gradient <- max(pmax(0, -gradient), pmax(0, gradient) * pmin(1, gap))
    c(coefficient = max(abs(problem$gr(par))), function_value = function_gradient)
  }
  append_fit <- function(method, fit) {
    par <- fit$par
    diagnostics <- assess(par)
    trace <<- rbind(trace, data.frame(method = method, convergence = fit$convergence,
      loss = problem$fn(par), coefficient_gradient = diagnostics[1L], function_gradient = diagnostics[2L]))
    fits[[length(fits) + 1L]] <<- list(par = par, optimizer = fit, diagnostics = diagnostics)
  }
  fit <- stats::optim(start, problem$fn, problem$gr, method = "BFGS",
                       control = list(maxit = max_iter, reltol = tolerance))
  append_fit("BFGS", fit)
  threshold <- 10 * sqrt(tolerance) * max(1, sqrt(abs(initial_loss)))
  valid <- function(fit) all(is.finite(fit$par)) &&
    all(is.finite(fit$diagnostics)) && max(fit$diagnostics) <= threshold
  if (!valid(fits[[1L]])) {
    fit <- stats::nlminb(fits[[1L]]$par, objective = problem$fn, gradient = problem$gr,
      control = list(iter.max = max_iter, eval.max = max_iter * 5L, rel.tol = tolerance,
                     x.tol = tolerance))
    append_fit("nlminb", fit)
  }
  for (restart in seq_len(8L)) {
    if (any(vapply(fits, valid, logical(1L)))) break
    best_index <- which.min(vapply(fits, function(fit) problem$fn(fit$par), numeric(1L)))
    candidate <- fits[[best_index]]$par
    gap <- exp(problem$sign * candidate)
    gradient <- 2 * as.numeric(S %*% problem$value(candidate) - b) / sqrt(diag(S))
    hidden_descent <- which(gradient < -threshold & gap < .1)
    if (!length(hidden_descent)) break
    # A finite restart revives a nearly flat link where increasing beta would
    # improve the objective. These are starting values, never coefficient
    # constraints: subsequent optimization is again over all real values.
    candidate[hidden_descent] <- problem$sign * log(.1)
    fit <- stats::optim(candidate, problem$fn, problem$gr, method = "BFGS",
                        control = list(maxit = max_iter, reltol = tolerance))
    append_fit(paste0("BFGS restart ", restart), fit)
  }
  eligible <- which(vapply(fits, valid, logical(1L)))
  if (!length(eligible)) {
    error <- structure(list(message = "Inverse-link sieve did not pass its convergence checks.",
      call = NULL, solver_trace = trace, parameters = fits[[which.min(trace$loss)]]$par,
      objective_matrix = S, objective_loading = b, objective_constant = constant, link = link), class = c("cmbridge_sieve_convergence_error", "error", "condition"))
    stop(error)
  }
  losses <- vapply(fits[eligible], function(fit) problem$fn(fit$par), numeric(1L))
  best <- fits[[eligible[which.min(losses)]]]
  list(values = problem$value(best$par), parameters = best$par, optimizer = best$optimizer,
       trace = trace, coefficient_gradient = unname(best$diagnostics[1L]),
       function_gradient = unname(best$diagnostics[2L]), tolerance = threshold,
       converged = TRUE, optimizer_reported_convergence = best$optimizer$convergence == 0L,
       accepted_by = "coefficient_and_function_gradient_checks", objective_scale = problem$scale,
       penalized_loss = problem$reference_loss + problem$fn(best$par) * problem$scale)
}
