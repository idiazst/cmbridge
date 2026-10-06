# Link parameters are unrestricted real coefficients. Bounds follow from the
# model parameterization; fitted functions and coefficients are never clipped.
.cm_link <- function(eta, link) {
  if (link == "inverse_logit") 1 + exp(-eta) else exp(eta)
}
.cm_link_derivative <- function(eta, link) {
  if (link == "inverse_logit") -exp(-eta) else exp(eta)
}

# Aggregate exact repeated target inputs, retaining every predictor. This is
# an algebraic reduction of the empirical moment, not a history restriction.
.linked_moment_design <- function(y, d, x, Hactive, Q, W) {
  active <- abs(d) > 0
  keys <- .cell_keys(x[active, , drop = FALSE])
  levels <- unique(keys); index <- match(keys, levels)
  first <- match(levels, keys)
  H <- Hactive[first, , drop = FALSE]
  assignment <- Matrix::sparseMatrix(i = seq_along(index), j = index,
    x = d[active], dims = c(length(index), length(levels)))
  K <- as.matrix(crossprod(Q[active, , drop = FALSE], assignment)) / length(y)
  cvec <- as.numeric(crossprod(Q, y) / length(y))
  list(H = H, K = K, cvec = cvec, W = W, index = index, active = active)
}

.linked_moment_problem <- function(design, link, penalty = 0) {
  H <- design$H; K <- design$K; cvec <- design$cvec; W <- design$W
  evaluate <- function(theta) {
    eta <- as.numeric(H %*% theta)
    value <- .cm_link(eta, link)
    if (any(!is.finite(value))) return(list(loss = Inf, gradient = rep(NA_real_, length(theta))))
    derivative <- .cm_link_derivative(eta, link)
    moment <- cvec - as.numeric(K %*% value)
    weighted <- as.numeric(W %*% moment)
    list(loss = as.numeric(crossprod(moment, weighted)) + penalty * sum(theta^2),
      gradient = 2 * (-as.numeric(crossprod(H, derivative * as.numeric(crossprod(K, weighted)))) + penalty * theta),
      moment = moment, values = value, derivative = derivative)
  }
  list(evaluate = evaluate, fn = function(theta) evaluate(theta)$loss,
       gr = function(theta) evaluate(theta)$gradient)
}

# Ridge sieve on a fixed basis with an unpenalized intercept. Scaling every
# nonconstant column is invertible; no predictors or interactions are removed.
.fit_linked_sieve <- function(y, d, x, Hactive, Q, W, xspec, zspec, ctrl) {
  scales <- sqrt(colMeans(Hactive^2)); scales[!is.finite(scales) | scales < 1e-8] <- 1
  scales[1L] <- 1
  penalty_weights <- rep(1, ncol(Hactive)); penalty_weights[1L] <- 0
  if (xspec$type == "cell_linear") {
    categories <- seq.int(xspec$p + 2L, ncol(Hactive))
    scales[categories] <- 1
    penalty_weights[categories] <- ctrl$cell_penalty_multiplier
  }
  design <- .linked_moment_design(y, d, x, sweep(Hactive, 2L, scales, "/"), Q, W)
  H <- design$H; K <- design$K; cvec <- design$cvec
  penalized <- penalty_weights
  evaluate <- function(theta) {
    eta <- as.numeric(H %*% theta)
    value <- .cm_link(eta, ctrl$link)
    if (any(!is.finite(value))) return(list(loss = Inf, gradient = rep(NA_real_, length(theta))))
    moment <- cvec - as.numeric(K %*% value)
    weighted <- as.numeric(W %*% moment)
    derivative <- .cm_link_derivative(eta, ctrl$link)
    list(loss = as.numeric(crossprod(moment, weighted)) + ctrl$lambda * sum(penalized * theta^2),
      gradient = -2 * as.numeric(crossprod(H, derivative * as.numeric(crossprod(K, weighted)))) +
        2 * ctrl$lambda * theta * penalized,
      values = value, moment = moment, derivative = derivative)
  }
  theta <- numeric(ncol(H))
  gap <- max(.1, mean(y) / mean(d) - if (ctrl$link == "inverse_logit") 1 else 0)
  theta[1L] <- if (ctrl$link == "inverse_logit") -log(gap) else log(gap)
  initial <- current <- evaluate(theta)
  threshold <- 10 * sqrt(ctrl$solver_tolerance) * max(1, sqrt(abs(initial$loss)))
  # Damped Gauss-Newton updates are unrestricted. The intercept is eliminated
  # by a Schur complement, retaining precision even at very large penalties.
  used <- 0L
  for (iter in seq_len(min(ctrl$solver_max_iter, 200L))) {
    if (max(abs(current$gradient)) <= threshold / 10) break
    J <- K %*% (current$derivative * H)
    curvature <- crossprod(J, W %*% J)
    curvature <- (curvature + t(curvature)) / 2
    curvature <- curvature + diag(ctrl$lambda * penalized + 1e-12, ncol(H))
    g <- current$gradient / 2
    if (ncol(H) == 1L) direction <- g / curvature[1L, 1L] else {
      block <- curvature[-1L, -1L, drop = FALSE]
      solved <- .safe_solve(block, cbind(g[-1L], curvature[-1L, 1L]))
      pivot <- curvature[1L, 1L] - sum(curvature[1L, -1L] * solved[, 2L])
      common <- (g[1L] - sum(curvature[1L, -1L] * solved[, 1L])) / max(pivot, 1e-12)
      direction <- c(common, solved[, 1L] - solved[, 2L] * common)
    }
    descent <- sum(g * direction)
    if (!is.finite(descent) || descent <= 0) direction <- g / pmax(diag(curvature), 1e-12)
    step <- 1
    for (attempt in 0:60) {
      proposal <- theta - step * direction
      candidate <- evaluate(proposal)
      if (is.finite(candidate$loss) && candidate$loss <= current$loss +
          8 * .Machine$double.eps * max(1, abs(current$loss))) break
      step <- step / 2
    }
    if (!is.finite(candidate$loss)) break
    theta <- proposal; current <- candidate; used <- iter
  }
  if (max(abs(current$gradient)) > threshold) {
    fn <- function(par) evaluate(par)$loss
    gr <- function(par) evaluate(par)$gradient
    parscale <- rep(1 / sqrt(max(ctrl$lambda, 1e-6)), length(theta)); parscale[1L] <- 1
    fit <- stats::optim(theta, fn, gr, method = "BFGS", control = list(
      maxit = ctrl$solver_max_iter, reltol = min(ctrl$solver_tolerance, 1e-14), parscale = parscale))
    if (fn(fit$par) <= current$loss) {theta <- fit$par; current <- evaluate(theta)}
  }
  if (!is.finite(current$loss) || max(abs(current$gradient)) > threshold)
    stop("Linked quadratic sieve did not pass its coefficient-gradient check.", call. = FALSE)
  coef <- theta / scales
  fitted <- rep(NA_real_, length(y)); fitted[design$active] <- current$values[design$index]
  list(coefficients = coef, link_coefficients = coef, target_spec = xspec, instrument_spec = zspec,
    target_design_cols = ncol(H), fitted = fitted, residual = y - d * ifelse(is.na(fitted), 0, fitted),
    moment_loss = as.numeric(crossprod(current$moment, W %*% current$moment)), tuning = ctrl,
    solver = list(converged = TRUE, coefficient_gradient = max(abs(current$gradient)),
      tolerance = threshold, iterations = used, initial_loss = initial$loss, final_loss = current$loss,
      method = "unrestricted damped Gauss-Newton", regularization = "ridge on scaled link coefficients"),
    predict_fun = .prediction_closure(function(newx) .cm_link(as.numeric(
      .eval_basis_spec(.as_matrix(newx), spec) %*% coefficient), link),
      list(spec = xspec, coefficient = coef, link = ctrl$link)),
    moment_weight = W,
    instrument_features = .prediction_closure(function(newz) .eval_basis_spec(newz, spec), list(spec = zspec)))
}

.fit_linked_landweber <- function(y, d, x, Hactive, Q, W, xspec, zspec, ctrl) {
  design <- .linked_moment_design(y, d, x, Hactive, Q, W)
  # An invertible change of coefficient coordinates preserves the entire basis
  # class while improving the conditioning of nonlinear Landweber iteration.
  transform <- diag(ncol(design$H))
  if (isTRUE(ctrl$precondition)) {
    covariance <- crossprod(Hactive) / nrow(Hactive)
    decomposition <- eigen((covariance + t(covariance)) / 2, symmetric = TRUE)
    floor <- max(decomposition$values) * 1e-3
    transform <- decomposition$vectors %*% diag(1 / sqrt(pmax(decomposition$values, floor)), ncol(Hactive)) %*% t(decomposition$vectors)
    design$H <- design$H %*% transform
  }
  problem <- .linked_moment_problem(design, ctrl$link)
  theta <- numeric(ncol(design$H))
  if (isTRUE(ctrl$init_intercept) && abs(mean(d)) > 1e-10) {
    value <- mean(y) / mean(d)
    gap <- max(.1, value - if (ctrl$link == "inverse_logit") 1 else 0)
    initial_coef <- numeric(ncol(design$H)); initial_coef[1L] <- if (ctrl$link == "inverse_logit") -log(gap) else log(gap)
    theta <- as.numeric(solve(transform, initial_coef))
  }
  n_iter <- if (is.null(ctrl$n_iter)) as.integer(ctrl$max_iter) else as.integer(ctrl$n_iter)
  if (n_iter < 1L) stop("n_iter must be positive.", call. = FALSE)
  if (!is.finite(ctrl$step_fraction) || ctrl$step_fraction <= 0 || ctrl$step_fraction > 1)
    stop("step_fraction must lie in (0,1].", call. = FALSE)
  initial <- current <- problem$evaluate(theta)
  jacobian <- design$K %*% (current$derivative * design$H)
  curvature <- crossprod(jacobian, W %*% jacobian)
  norm <- max(eigen((curvature + t(curvature)) / 2, symmetric = TRUE, only.values = TRUE)$values)
  if (!is.finite(norm) || norm <= 0)
    stop("Conditional-moment operator has zero numerical norm.", call. = FALSE)
  maximum_step <- ctrl$step_fraction / norm
  # A fixed initial step with backtracking is nonlinear Landweber iteration.
  # Re-inverting curvature at every iterate can explode the step when the
  # inverse-expit derivative approaches zero in a small training sample.
  # Backtracking may reduce this step, but never increases it. Coefficients
  # remain unrestricted and iteration stopping remains the regularization.
  used <- 0L; last_delta <- Inf; step <- maximum_step; reductions <- 0L
  for (iter in seq_len(n_iter)) {
    # Half the squared-loss gradient matches the original Landweber convention.
    gradient <- current$gradient / 2
    if (max(abs(gradient)) <= ctrl$tol) break
    accepted <- FALSE
    for (attempt in 0:60) {
      proposal <- theta - step * gradient
      candidate <- problem$evaluate(proposal)
      if (is.finite(candidate$loss) && candidate$loss <=
          current$loss - 1e-4 * step * sum(gradient^2)) {accepted <- TRUE; break}
      step <- step / 2
      reductions <- reductions + 1L
    }
    if (!accepted) stop("Linked Landweber could not find a decreasing unrestricted update.", call. = FALSE)
    last_delta <- max(abs(proposal - theta))
    theta <- proposal; current <- candidate; used <- iter
    if (last_delta < ctrl$tol) break
  }
  pred_train <- rep(NA_real_, length(y))
  pred_train[design$active] <- current$values[design$index]
  residual <- y - d * ifelse(is.na(pred_train), 0, pred_train)
  moment <- as.numeric(crossprod(Q, residual) / length(y))
  coefficient <- as.numeric(transform %*% theta)
  list(coefficients = coefficient, link_coefficients = coefficient, target_spec = xspec, instrument_spec = zspec,
    fitted = pred_train, residual = residual, moment_loss = as.numeric(crossprod(moment, W %*% moment)),
    tuning = c(ctrl, list(step = step, maximum_step = maximum_step,
      iterations_used = used, final_delta = last_delta)),
    solver = list(method = "nonlinear Landweber", initial_loss = initial$loss, final_loss = current$loss,
      iterations = used, maximum_step = maximum_step, final_step = step,
      backtracking_reductions = reductions, coefficient_gradient = max(abs(current$gradient)),
      stopped_by_tolerance = last_delta < ctrl$tol || max(abs(current$gradient / 2)) <= ctrl$tol,
      regularization = "iteration stopping"),
    predict_fun = .prediction_closure(function(newx) {
      eta <- as.numeric(.eval_basis_spec(.as_matrix(newx), spec) %*% coefficient)
      .cm_link(eta, link)
    }, list(spec = xspec, coefficient = coefficient, link = ctrl$link)),
    moment_weight = W,
    instrument_features = .prediction_closure(function(newz) .eval_basis_spec(newz, spec), list(spec = zspec)))
}

.fit_linked_pmmr <- function(y, d, x, Phi_active, Psi, ctrl, lambda) {
  W <- diag(ncol(Psi))
  design <- .linked_moment_design(y, d, x, Phi_active, Psi, W)
  problem <- .linked_moment_problem(design, ctrl$link, penalty = lambda)
  start <- numeric(ncol(Phi_active))
  initial <- problem$fn(start)
  fit <- stats::optim(start, problem$fn, problem$gr, method = "BFGS",
    control = list(maxit = ctrl$solver_max_iter, reltol = ctrl$solver_tolerance))
  trace <- data.frame(method = "BFGS", stop_code = fit$convergence, loss = problem$fn(fit$par),
    coefficient_gradient = max(abs(problem$gr(fit$par))))
  threshold <- 10 * sqrt(ctrl$solver_tolerance) * max(1, sqrt(abs(initial)))
  if (trace$coefficient_gradient > threshold) {
    refined <- stats::nlminb(fit$par, objective = problem$fn, gradient = problem$gr,
      control = list(iter.max = ctrl$solver_max_iter, eval.max = ctrl$solver_max_iter * 5L,
        rel.tol = ctrl$solver_tolerance, x.tol = ctrl$solver_tolerance))
    trace <- rbind(trace, data.frame(method = "nlminb", stop_code = refined$convergence,
      loss = problem$fn(refined$par), coefficient_gradient = max(abs(problem$gr(refined$par)))))
    if (tail(trace$loss, 1L) <= trace$loss[1L]) fit <- refined
  }
  final <- problem$evaluate(fit$par)
  if (any(!is.finite(fit$par)) || !is.finite(final$loss) || max(abs(final$gradient)) > threshold)
    stop("Linked PMMR did not pass its coefficient-gradient check.", call. = FALSE)
  pred_train <- rep(NA_real_, length(y))
  pred_train[design$active] <- final$values[design$index]
  list(coefficients = fit$par, link_coefficients = fit$par, fitted = pred_train,
    residual = y - d * ifelse(is.na(pred_train), 0, pred_train), moment_loss = sum(final$moment^2),
    solver = list(method = "unrestricted nonlinear PMMR", initial_loss = initial, final_loss = final$loss,
      penalized_loss = final$loss, coefficient_gradient = max(abs(final$gradient)), tolerance = threshold,
      optimizer_stop_code = fit$convergence, trace = trace, converged = TRUE),
    regularization = "RKHS penalty on the linear predictor")
}
