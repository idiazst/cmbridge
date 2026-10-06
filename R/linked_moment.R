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

.fit_linked_landweber <- function(y, d, x, Hactive, Q, W, xspec, zspec, ctrl) {
  design <- .linked_moment_design(y, d, x, Hactive, Q, W)
  problem <- .linked_moment_problem(design, ctrl$link)
  theta <- numeric(ncol(design$H))
  if (isTRUE(ctrl$init_intercept) && abs(mean(d)) > 1e-10) {
    value <- mean(y) / mean(d)
    gap <- max(.1, value - if (ctrl$link == "inverse_logit") 1 else 0)
    theta[1L] <- if (ctrl$link == "inverse_logit") -log(gap) else log(gap)
  }
  n_iter <- if (is.null(ctrl$n_iter)) as.integer(ctrl$max_iter) else as.integer(ctrl$n_iter)
  if (n_iter < 1L) stop("n_iter must be positive.", call. = FALSE)
  if (!is.finite(ctrl$step_fraction) || ctrl$step_fraction <= 0 || ctrl$step_fraction > 1)
    stop("step_fraction must lie in (0,1].", call. = FALSE)
  initial <- current <- problem$evaluate(theta)
  used <- 0L; last_delta <- Inf; step <- NA_real_; reductions <- 0L
  for (iter in seq_len(n_iter)) {
    jacobian <- design$K %*% (current$derivative * design$H)
    curvature <- crossprod(jacobian, W %*% jacobian)
    norm <- max(eigen((curvature + t(curvature)) / 2, symmetric = TRUE, only.values = TRUE)$values)
    if (!is.finite(norm) || norm <= 0) break
    # Half the squared-loss gradient matches the original Landweber convention.
    gradient <- current$gradient / 2
    if (max(abs(gradient)) <= ctrl$tol) break
    step <- ctrl$step_fraction / norm
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
  list(coefficients = theta, link_coefficients = theta, target_spec = xspec, instrument_spec = zspec,
    fitted = pred_train, residual = residual, moment_loss = as.numeric(crossprod(moment, W %*% moment)),
    tuning = c(ctrl, list(step = step, iterations_used = used, final_delta = last_delta)),
    solver = list(method = "nonlinear Landweber", initial_loss = initial$loss, final_loss = current$loss,
      iterations = used, backtracking_reductions = reductions, coefficient_gradient = max(abs(current$gradient)),
      stopped_by_tolerance = last_delta < ctrl$tol || max(abs(current$gradient / 2)) <= ctrl$tol,
      regularization = "iteration stopping"),
    predict_fun = .prediction_closure(function(newx) {
      eta <- as.numeric(.eval_basis_spec(.as_matrix(newx), spec) %*% coefficient)
      .cm_link(eta, link)
    }, list(spec = xspec, coefficient = theta, link = ctrl$link)),
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
