.nystrom_spec <- function(centers, bandwidth, ridge = 1e-10, tol = 1e-10) {
  centers <- .as_matrix(centers)
  Kcc <- rbf_kernel(centers, centers, bandwidth)
  Kcc <- (Kcc + t(Kcc)) / 2
  ee <- eigen(Kcc, symmetric = TRUE)
  vmax <- max(ee$values)
  keep <- ee$values > max(tol * vmax, ridge)
  if (!any(keep)) stop("Nystrom kernel matrix has no numerically positive eigenvalues.", call. = FALSE)
  transform <- ee$vectors[, keep, drop = FALSE] %*%
    diag(1 / sqrt(ee$values[keep] + ridge), nrow = sum(keep))
  list(centers = centers, bandwidth = bandwidth, transform = transform)
}

.nystrom_apply <- function(x, spec) {
  rbf_kernel(.as_matrix(x), spec$centers, spec$bandwidth) %*% spec$transform
}

.pmmr_core <- function(y, d, x, z, ctrl, lambda) {
  active <- abs(d) > 0
  if (!any(active)) stop("diagonal is zero for every observation.", call. = FALSE)
  if (anyNA(x[active, , drop = FALSE])) stop("target may be missing only where diagonal is zero.", call. = FALSE)
  if (anyNA(z)) stop("instrument cannot be missing for PMMR.", call. = FALSE)

  xs <- x[active, , drop = FALSE]
  zs <- z
  xscale <- if (isTRUE(ctrl$scale)) .scale_fit(xs) else list(center = rep(0, ncol(xs)), scale = rep(1, ncol(xs)))
  zscale <- if (isTRUE(ctrl$scale)) .scale_fit(zs) else list(center = rep(0, ncol(zs)), scale = rep(1, ncol(zs)))
  xs_sc <- .scale_apply(xs, xscale)
  z_sc <- .scale_apply(zs, zscale)

  xc <- ctrl$target_centers_matrix
  if (is.null(xc)) xc <- .select_centers(xs_sc, ctrl$target_centers, ctrl$seed)
  else xc <- .scale_apply(.as_matrix(xc), xscale)
  zc <- ctrl$instrument_centers_matrix
  if (is.null(zc)) zc <- .select_centers(z_sc, ctrl$instrument_centers, ctrl$seed + 1L)
  else zc <- .scale_apply(.as_matrix(zc), zscale)

  bwx <- ctrl$target_bandwidth
  if (is.null(bwx)) bwx <- .median_distance(xs_sc, seed = ctrl$seed)
  bwz <- ctrl$instrument_bandwidth
  if (is.null(bwz)) bwz <- .median_distance(z_sc, seed = ctrl$seed + 1L)
  if (bwx <= 0 || bwz <= 0) stop("kernel bandwidths must be positive.", call. = FALSE)

  xspec <- .nystrom_spec(xc, bwx, ridge = ctrl$nystrom_ridge, tol = ctrl$nystrom_tol)
  zspec <- .nystrom_spec(zc, bwz, ridge = ctrl$nystrom_ridge, tol = ctrl$nystrom_tol)

  ## PMMR uses the V-statistic maximum-moment loss
  ##   n^{-2} r^T K_Z r + lambda ||h||_{H_X}^2.
  ## Under the Nyström approximations K_Z ~= Psi Psi^T and
  ## K_X ~= Phi Phi^T, this becomes
  ##   || n^{-1} Psi^T (y - D Phi theta) ||^2 + lambda ||theta||^2.
  ## The whitening below is only the Nyström RKHS normalization K_CC^{-1/2};
  ## unlike the previous implementation, there is no inverse covariance
  ## weighting of critic features.
  Phi_active <- .nystrom_apply(xs_sc, xspec)
  Phi <- matrix(0, nrow(x), ncol(Phi_active))
  Phi[active, ] <- Phi_active
  Psi <- .nystrom_apply(z_sc, zspec)
  n <- length(y)
  linked <- NULL
  if (ctrl$link != "identity") {
    linked <- .fit_linked_pmmr(y, d, x, Phi_active, Psi, ctrl, lambda)
    coef <- linked$coefficients
    pred_train <- linked$fitted
    residual <- linked$residual
  } else {
    A <- crossprod(Psi, d * Phi) / n
    b <- as.numeric(crossprod(Psi, y) / n)
    lhs <- crossprod(A) + lambda * diag(ncol(Phi))
    rhs <- crossprod(A, b)
    coef <- as.numeric(.safe_solve(lhs, rhs, ridge = ctrl$solve_ridge))
    pred_active <- as.numeric(Phi_active %*% coef)
    pred_train <- rep(NA_real_, n)
    pred_train[active] <- pred_active
    residual <- y - d * ifelse(is.na(pred_train), 0, pred_train)
  }
  moment <- as.numeric(crossprod(Psi, residual) / n)
  predict_fun <- .prediction_closure(function(newx) {
    eta <- as.numeric(.nystrom_apply(.scale_apply(.as_matrix(newx), scale), spec) %*% coefficient)
    if (link == "identity") eta else .cm_link(eta, link)
  }, list(scale = xscale, spec = xspec, coefficient = coef, link = ctrl$link))
  instrument_fun <- .prediction_closure(function(newz)
    .nystrom_apply(.scale_apply(.as_matrix(newz), scale), spec),
    list(scale = zscale, spec = zspec))

  list(
    coefficients = coef,
    link_coefficients = if (is.null(linked)) NULL else coef,
    solver = if (is.null(linked)) NULL else linked$solver,
    regularization = if (is.null(linked)) "RKHS penalty on function" else linked$regularization,
    target_centers = xc,
    instrument_centers = zc,
    target_scale = xscale,
    instrument_scale = zscale,
    target_bandwidth = bwx,
    instrument_bandwidth = bwz,
    target_nystrom = xspec,
    instrument_nystrom = zspec,
    fitted = pred_train,
    residual = residual,
    moment_loss = sum(moment^2),
    tuning = utils::modifyList(ctrl, list(lambda = lambda)),
    predict_fun = predict_fun,
    moment_weight = diag(ncol(Psi)),
    instrument_features = instrument_fun
  )
}

.fit_pmmr <- function(y, d, x, z, control) {
  ctrl <- .merge_control(list(
    target_centers = 40L, instrument_centers = 80L,
    target_bandwidth = NULL, instrument_bandwidth = NULL,
    target_centers_matrix = NULL, instrument_centers_matrix = NULL,
    lambda = 1e-4,
    nystrom_ridge = 1e-10, nystrom_tol = 1e-10,
    solve_ridge = 1e-12,
    scale = TRUE, seed = 1L, link = "identity",
    solver_tolerance = 1e-10, solver_max_iter = 10000L
  ), control)

  ctrl$link <- match.arg(ctrl$link, c("identity", "log", "inverse_logit"))
  if (!is.finite(ctrl$solver_tolerance) || ctrl$solver_tolerance <= 0 ||
      !is.finite(ctrl$solver_max_iter) || ctrl$solver_max_iter < 1)
    stop("Invalid linked-PMMR solver controls.", call. = FALSE)
  if (length(ctrl$lambda) != 1L || !is.finite(ctrl$lambda) || ctrl$lambda <= 0) {
    stop("lambda must be a single positive finite number.", call. = FALSE)
  }
  .pmmr_core(y, d, x, z, ctrl, lambda = ctrl$lambda)
}
