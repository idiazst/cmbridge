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

.fit_basis_spec <- function(x, type = c("bs", "poly", "cell"), degree = 3L, df = 6L) {
  type <- match.arg(type)
  x <- .as_matrix(x)
  if (anyNA(x)) stop("basis fitting data cannot contain missing values.", call. = FALSE)
  if (type == "cell") {
    return(list(type = type, levels = sort(unique(.cell_keys(x))), p = ncol(x)))
  }
  if (type == "poly") {
    return(list(type = type, degree = as.integer(degree), p = ncol(x)))
  }
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
  if (spec$type == "cell") {
    index <- match(.cell_keys(x), spec$levels)
    out <- matrix(0, nrow(x), length(spec$levels) + 1L)
    out[, 1L] <- 1
    seen <- which(!is.na(index))
    out[cbind(seen, 1L + index[seen])] <- 1
    return(out)
  }
  if (spec$type == "poly") return(poly_basis(x, degree = spec$degree, intercept = TRUE))
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
