#' Polynomial basis
#'
#' @param x Numeric vector or matrix.
#' @param degree Maximum marginal polynomial degree.
#' @param intercept Include an intercept column.
#' @return Numeric design matrix.
#' @export
poly_basis <- function(x, degree = 3L, intercept = TRUE) {
  x <- .as_matrix(x)
  degree <- as.integer(degree)
  if (degree < 1L) stop("degree must be at least 1.", call. = FALSE)
  out <- if (isTRUE(intercept)) matrix(1, nrow(x), 1L) else NULL
  for (j in seq_len(ncol(x))) {
    for (k in seq_len(degree)) out <- cbind(out, x[, j]^k)
  }
  colnames(out) <- NULL
  out
}

.fit_basis_spec <- function(x, type = c("bs", "poly", "cell", "quadratic", "cell_linear"), degree = 3L, df = 6L) {
  type <- match.arg(type)
  x <- .as_matrix(x)
  if (anyNA(x)) stop("basis fitting data cannot contain missing values.", call. = FALSE)
  if (type %in% c("cell", "cell_linear")) {
    return(list(type = type, levels = sort(unique(.cell_keys(x))), p = ncol(x)))
  }
  if (type == "poly") {
    return(list(type = type, degree = as.integer(degree), p = ncol(x)))
  }
  if (type == "quadratic") return(list(type = type, degree = 2L, p = ncol(x)))
  specs <- vector("list", ncol(x))
  for (j in seq_len(ncol(x))) {
    if (length(unique(x[, j])) == 1L) {
      # A constant column cannot identify a spline effect. Keep its learned
      # contribution constant for new values instead of evaluating splines
      # with identical boundary knots, which can return NaNs.
      specs[[j]] <- list(constant = TRUE)
      next
    }
    b <- splines::bs(x[, j], df = df, degree = degree, intercept = FALSE)
    knots <- attr(b, "knots"); boundary <- attr(b, "Boundary.knots")
    # bs() leaves all-equal boundary knots in place. Its extrapolation then
    # evaluates derivatives at a zero-width pivot and can return NaNs.
    # Use the same one-eighth inward spacing as bs()'s ordinary adjustment.
    if (length(knots) && (all(knots == boundary[1L]) || all(knots == boundary[2L]))) {
      knots[] <- if (all(knots == boundary[1L])) boundary[1L] + diff(boundary) / 8 else
        boundary[2L] - diff(boundary) / 8
      b <- splines::bs(x[, j], knots = knots, Boundary.knots = boundary,
                       degree = degree, intercept = FALSE)
    }
    specs[[j]] <- list(
      knots = attr(b, "knots"),
      Boundary.knots = attr(b, "Boundary.knots"),
      degree = attr(b, "degree")
    )
  }
  list(type = type, specs = specs, p = ncol(x))
}

.eval_basis_spec <- function(x, spec) {
  x <- .as_matrix(x)
  if (ncol(x) != spec$p) stop("new data have the wrong number of columns.", call. = FALSE)
  if (spec$type %in% c("cell", "cell_linear")) {
    index <- match(.cell_keys(x), spec$levels)
    out <- matrix(0, nrow(x), length(spec$levels) + 1L)
    out[, 1L] <- 1
    seen <- which(!is.na(index))
    out[cbind(seen, 1L + index[seen])] <- 1
    if (spec$type == "cell_linear") return(cbind(poly_basis(x, degree = 1L), out[, -1L, drop = FALSE]))
    return(out)
  }
  if (spec$type == "poly") return(poly_basis(x, degree = spec$degree, intercept = TRUE))
  if (spec$type == "quadratic") {
    # A fixed polynomial basis: intercept, every main effect, every square,
    # and every pairwise interaction. No levels or terms are learned from x.
    out <- poly_basis(x, degree = 2L, intercept = TRUE)
    if (spec$p > 1L) for (j in seq_len(spec$p - 1L)) for (k in (j + 1L):spec$p)
      out <- cbind(out, x[, j] * x[, k])
    return(out)
  }
  if (!nrow(x)) {
    # A bridge validation sample may have no measured outcomes. Its spline
    # prediction design is empty, but must retain the fitted column count.
    # splines::bs() does not accept a zero-length vector with fixed knots.
    width <- 1L + sum(vapply(spec$specs, function(sj) {
      if (isTRUE(sj$constant)) 0L else length(sj$knots) + as.integer(sj$degree)
    }, integer(1)))
    return(matrix(numeric(), nrow = 0L, ncol = width))
  }
  out <- matrix(1, nrow(x), 1L)
  for (j in seq_len(ncol(x))) {
    sj <- spec$specs[[j]]
    if (isTRUE(sj$constant)) next
    value <- x[, j]
    if (identical(spec$extrapolation, "constant")) {
      # Constant continuation is part of the basis definition. It never
      # clips a fitted function or its link coefficients.
      value <- pmin(sj$Boundary.knots[2L], pmax(sj$Boundary.knots[1L], value))
    }
    bj <- splines::bs(
      value, knots = sj$knots, Boundary.knots = sj$Boundary.knots,
      degree = sj$degree, intercept = FALSE
    )
    out <- cbind(out, bj)
  }
  out
}

# Remove the redundant intercept from an indicator critic without changing its
# ridge-weighted moment loss. F=(1,I) maps the nonredundant indicators to the
# original design; F(C+ridge I)^-1 F' is evaluated by Sherman-Morrison.
.moment_basis_spec <- function(z, spec, ridge) {
  Q <- .eval_basis_spec(z, spec)
  if (spec$type == "cell") {
    Q <- Q[, -1L, drop = FALSE]
    counts <- colMeans(Q)
    inverse <- 1 / (counts + ridge)
    alpha <- ridge / (length(counts) + 1)
    W <- diag(inverse, length(counts)) +
      alpha / (1 - alpha * sum(inverse)) * tcrossprod(inverse)
  } else W <- .safe_inverse(crossprod(Q) / nrow(Q) + ridge * diag(ncol(Q)), ridge)
  features <- .prediction_closure(function(newz) {
    out <- .eval_basis_spec(.as_matrix(newz), spec)
    if (spec$type == "cell") out[, -1L, drop = FALSE] else out
  }, list(spec = spec))
  list(Q = Q, W = W, features = features)
}

# The common cell value is not penalized. Eliminate it separately so a large
# ridge cannot erase its small curvature through floating-point cancellation.
.centered_ridge_solve <- function(S, b, lambda) {
  p <- nrow(S)
  if (p == 1L) return(as.numeric(.safe_solve(S, b)))
  basis <- qr.Q(qr(cbind(rep(1 / sqrt(p), p), diag(p)[, -1L, drop = FALSE])))
  small <- crossprod(basis, S %*% basis)
  rhs <- as.numeric(crossprod(basis, b))
  block <- small[-1L, -1L, drop = FALSE] + diag(lambda, p - 1L)
  solved <- .safe_solve(block, cbind(rhs[-1L], small[-1L, 1L]))
  pivot <- small[1L, 1L] - sum(small[1L, -1L] * solved[, 2L])
  common <- (rhs[1L] - sum(small[1L, -1L] * solved[, 1L])) / pivot
  as.numeric(basis %*% c(common, solved[, 1L] - solved[, 2L] * common))
}

#' Gaussian radial-basis kernel
#'
#' @param x Numeric vector or matrix.
#' @param centers Numeric vector or matrix of kernel centers.
#' @param bandwidth Positive Gaussian bandwidth.
#' @return Kernel design matrix.
#' @export
rbf_kernel <- function(x, centers, bandwidth) {
  x <- .as_matrix(x)
  centers <- .as_matrix(centers)
  if (ncol(x) != ncol(centers)) stop("x and centers must have the same number of columns.", call. = FALSE)
  if (!is.finite(bandwidth) || bandwidth <= 0) stop("bandwidth must be positive.", call. = FALSE)
  x2 <- rowSums(x^2)
  c2 <- rowSums(centers^2)
  d2 <- outer(x2, c2, "+") - 2 * tcrossprod(x, centers)
  d2[d2 < 0] <- 0
  exp(-d2 / (2 * bandwidth^2))
}
