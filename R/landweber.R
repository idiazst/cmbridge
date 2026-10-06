.fit_landweber <- function(y, d, x, z, control) {
  ctrl <- .merge_control(list(
    target_basis = "bs", instrument_basis = "bs",
    target_degree = 3L, instrument_degree = 3L,
    target_df = 6L, instrument_df = 9L,
    max_iter = 2000L, n_iter = NULL,
    step_fraction = 0.95, weight_ridge = 1e-8,
    init_intercept = TRUE, tol = 1e-10, link = "identity", target_extrapolation = NULL,
    precondition = FALSE
  ), control)

  ctrl$link <- match.arg(ctrl$link, c("identity", "log", "inverse_logit"))
  if (is.null(ctrl$target_extrapolation))
    ctrl$target_extrapolation <- if (ctrl$link == "identity") "polynomial" else "constant"
  ctrl$target_extrapolation <- match.arg(ctrl$target_extrapolation, c("constant", "polynomial"))
  active <- abs(d) > 0
  if (!any(active)) stop("diagonal is zero for every observation.", call. = FALSE)
  if (anyNA(x[active, , drop = FALSE])) stop("target may be missing only where diagonal is zero.", call. = FALSE)
  if (anyNA(z)) stop("instrument cannot be missing for Landweber estimation.", call. = FALSE)

  xspec <- .fit_basis_spec(x[active, , drop = FALSE], ctrl$target_basis,
                           ctrl$target_degree, ctrl$target_df)
  if (xspec$type == "bs") xspec$extrapolation <- ctrl$target_extrapolation
  zspec <- .fit_basis_spec(z, ctrl$instrument_basis,
                           ctrl$instrument_degree, ctrl$instrument_df)
  Hactive <- .eval_basis_spec(x[active, , drop = FALSE], xspec)
  H <- matrix(0, nrow(x), ncol(Hactive)); H[active, ] <- Hactive
  critic <- .moment_basis_spec(z, zspec, ctrl$weight_ridge)
  Q <- critic$Q
  if (ctrl$link != "identity") {
    out <- .fit_linked_landweber(y, d, x, Hactive, Q, critic$W, xspec, zspec, ctrl)
    out$instrument_features <- critic$features
    return(out)
  }
  n <- length(y)

  A <- crossprod(Q, d * H) / n
  cvec <- crossprod(Q, y) / n
  W <- critic$W
  S <- crossprod(A, W %*% A)
  b <- as.numeric(crossprod(A, W %*% cvec))
  eigmax <- max(eigen((S + t(S)) / 2, symmetric = TRUE, only.values = TRUE)$values)
  if (!is.finite(eigmax) || eigmax <= 0) stop("conditional-moment operator has zero numerical norm.", call. = FALSE)
  step <- ctrl$step_fraction / eigmax

  theta <- rep(0, ncol(H))
  if (isTRUE(ctrl$init_intercept) && abs(mean(d)) > 1e-10) theta[1L] <- mean(y) / mean(d)
  n_iter <- if (is.null(ctrl$n_iter)) as.integer(ctrl$max_iter) else as.integer(ctrl$n_iter)
  if (n_iter < 1L) stop("n_iter must be positive.", call. = FALSE)

  last_delta <- Inf
  used <- 0L
  for (iter in seq_len(n_iter)) {
    grad <- b - as.numeric(S %*% theta)
    next_theta <- theta + step * grad
    last_delta <- max(abs(next_theta - theta))
    theta <- next_theta
    used <- iter
    if (last_delta < ctrl$tol) break
  }

  pred_train <- as.numeric(H %*% theta)
  pred_train[!active] <- NA_real_
  residual <- y - d * ifelse(is.na(pred_train), 0, pred_train)
  moment <- as.numeric(crossprod(Q, residual) / n)

  list(
    coefficients = theta,
    target_spec = xspec,
    instrument_spec = zspec,
    fitted = pred_train,
    residual = residual,
    moment_loss = as.numeric(crossprod(moment, W %*% moment)),
    tuning = c(ctrl, list(step = step, iterations_used = used, final_delta = last_delta)),
    predict_fun = .prediction_closure(function(newx)
      as.numeric(.eval_basis_spec(.as_matrix(newx), spec) %*% coefficient),
      list(spec = xspec, coefficient = theta)),
    moment_weight = W,
    instrument_features = critic$features
  )
}
