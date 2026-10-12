#' Fit a conditional-moment function
#'
#' Fits a function h satisfying E[Y - D h(X) | Z] = 0 using one of four
#' standardized learners.
#'
#' @param response Numeric response Y.
#' @param diagonal Numeric diagonal/multiplicative variable D.
#' @param target Numeric vector or matrix X, the argument of h. Values may be
#'   missing only on rows where diagonal is zero.
#' @param instrument Numeric vector or matrix Z defining the conditional moment.
#' @param method One of "sieve_md", "landweber", "pmmr", or "saturated_l1".
#' @param control Method-specific control list. The cell sieve, Landweber and PMMR support
#'   positive links with unrestricted real coefficients. Optional lower and upper
#'   bounds must be infinite; post-fit clipping is not supported.
#' @return An object of class cmbridge_fit.
#' @export
fit_cm <- function(response, diagonal = NULL, target, instrument,
                   method = c("sieve_md", "landweber", "pmmr", "saturated_l1"), control = list()) {
  method <- match.arg(method)
  y <- as.numeric(response)
  d <- if (is.null(diagonal)) rep(1, length(y)) else as.numeric(diagonal)
  x <- .as_matrix(target, "target")
  z <- .as_matrix(instrument, "instrument")
  n <- .check_lengths(y, d, x, z)
  if (any(!is.finite(y))) stop("response must be finite.", call. = FALSE)
  if (any(!is.finite(d))) stop("diagonal must be finite.", call. = FALSE)
  if (anyNA(z) || any(!is.finite(z))) stop("instrument must be finite and nonmissing.", call. = FALSE)
  active <- abs(d) > 0
  if (any(!is.finite(x[active, , drop = FALSE]))) stop("target must be finite wherever diagonal is nonzero.", call. = FALSE)

  lower <- if (is.null(control$lower)) -Inf else control$lower
  upper <- if (is.null(control$upper)) Inf else control$upper
  if (length(lower) != 1L || length(upper) != 1L ||
      is.na(lower) || is.na(upper) || lower > upper) {
    stop("lower and upper must be ordered scalar prediction bounds.", call. = FALSE)
  }
  if (is.finite(lower) || is.finite(upper)) {
    stop("Post-fit prediction bounds are not supported. Use a supported model parameterization instead of post-fit bounds.", call. = FALSE)
  }
  core <- switch(method,
    sieve_md = .fit_sieve_md(y, d, x, z, control),
    landweber = .fit_landweber(y, d, x, z, control),
    pmmr = .fit_pmmr(y, d, x, z, control),
    saturated_l1 = .fit_saturated_l1(y, d, x, z, control)
  )

  out <- c(core, list(
    method = method,
    n = n,
    response = y,
    diagonal = d,
    target_dim = ncol(x),
    instrument_dim = ncol(z),
    call = match.call(),
    type = "generic"
  ))
  class(out) <- "cmbridge_fit"
  out
}

#' Fit an inverse-probability bridge
#'
#' Fits beta in E[M beta(V) - 1 | B] = 0.
#'
#' @param B Shadow-side variables.
#' @param V Bridge-side variables. May be missing when M=0.
#' @param M Response/visit indicator.
#' @param method Estimation method.
#' @param control Method-specific control list.
#' @return A cmbridge_fit object.
#' @export
fit_bridge <- function(B, V, M, method = c("sieve_md", "landweber", "pmmr", "saturated_l1"), control = list()) {
  M <- as.numeric(M)
  if (anyNA(M) || any(M < 0) || any(M > 1)) stop("M must lie in [0,1].", call. = FALSE)
  fit <- fit_cm(rep(1, length(M)), M, V, B, method = match.arg(method), control = control)
  fit$type <- "bridge"
  fit
}

#' Fit an adjoint representer
#'
#' Fits lambda in E[M {lambda(B) - phi(V)} | V] = 0. Since V is only
#' required on complete cases, estimation is carried out among M>0, where this
#' is the equivalent conditional moment E[lambda(B)-phi(V)|V,M>0]=0.
#'
#' @param B Variables indexing lambda.
#' @param V Variables conditioning the adjoint equation.
#' @param M Response/visit indicator.
#' @param phi Numeric target-specific loading evaluated on complete cases.
#' @param method Estimation method.
#' @param control Method-specific control list.
#' @return A cmbridge_fit object.
#' @export
fit_adjoint <- function(B, V, M, phi, method = c("sieve_md", "landweber", "pmmr", "saturated_l1"), control = list()) {
  B <- .as_matrix(B, "B")
  V <- .as_matrix(V, "V")
  M <- as.numeric(M)
  phi <- as.numeric(phi)
  if (length(M) != nrow(B) || nrow(V) != nrow(B) || length(phi) != nrow(B)) {
    stop("B, V, M, and phi must have the same number of observations.", call. = FALSE)
  }
  keep <- is.finite(M) & M > 0
  if (!any(keep)) stop("no complete cases are available for adjoint estimation.", call. = FALSE)
  if (any(!is.finite(phi[keep]))) stop("phi must be finite on complete cases.", call. = FALSE)
  if (any(!is.finite(B[keep, , drop = FALSE])) || any(!is.finite(V[keep, , drop = FALSE]))) {
    stop("B and V must be finite on complete cases.", call. = FALSE)
  }
  fit <- fit_cm(phi[keep], rep(1, sum(keep)), B[keep, , drop = FALSE],
                V[keep, , drop = FALSE], method = match.arg(method), control = control)
  fit$type <- "adjoint"
  fit$complete_case_index <- which(keep)
  fit
}

#' @export
predict.cmbridge_fit <- function(object, newdata, ...) {
  object$predict_fun(newdata)
}

#' @export
fitted.cmbridge_fit <- function(object, ...) object$fitted

#' @export
print.cmbridge_fit <- function(x, ...) {
  cat("Conditional-moment fit\n")
  cat("  type:   ", x$type, "\n", sep = "")
  cat("  method: ", x$method, "\n", sep = "")
  cat("  n:      ", x$n, "\n", sep = "")
  cat("  training moment loss: ", format(x$moment_loss, digits = 5), "\n", sep = "")
  invisible(x)
}

#' Evaluate a finite-feature conditional-moment loss
#'
#' @param object Fitted cmbridge_fit object.
#' @param response Optional response vector; defaults to training response.
#' @param diagonal Optional diagonal vector; defaults to training diagonal.
#' @param target Target values for prediction. Required for new data.
#' @param instrument Instrument values. Required for new data.
#' @return Squared Euclidean norm of the empirical feature moments.
#' @export
moment_loss <- function(object, response = NULL, diagonal = NULL,
                        target = NULL, instrument = NULL) {
  if (!inherits(object, "cmbridge_fit")) stop("object must be a cmbridge_fit.", call. = FALSE)
  if (is.null(response) && is.null(diagonal) && is.null(target) && is.null(instrument)) return(object$moment_loss)
  if (is.null(response) || is.null(diagonal) || is.null(target) || is.null(instrument)) {
    stop("response, diagonal, target, and instrument must all be supplied for new data.", call. = FALSE)
  }
  y <- as.numeric(response); d <- as.numeric(diagonal)
  x <- .as_matrix(target); z <- .as_matrix(instrument)
  .check_lengths(y, d, x, z)
  active <- abs(d) > 0
  pred <- rep(0, length(y))
  if (any(active)) pred[active] <- predict(object, x[active, , drop = FALSE])
  r <- y - d * pred
  if (inherits(object, "cmbridge_ensemble")) {
    if (any(!is.finite(z)) || any(!is.finite(r))) {
      stop("new-data moment loss requires finite instruments and residuals.", call. = FALSE)
    }
    return(.ensemble_gram(matrix(r, ncol = 1L), z, object$scoring_kernel)[1, 1])
  }
  Q <- object$instrument_features(z)
  m <- as.numeric(crossprod(Q, r) / length(y))
  W <- object$moment_weight
  if (is.null(W) || nrow(W) != length(m)) return(sum(m^2))
  as.numeric(crossprod(m, W %*% m))
}
