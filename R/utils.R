.prediction_closure <- function(fun, bindings) {
  # Retain prediction parameters without retaining the fitting frame's arrays.
  environment(fun) <- list2env(bindings, parent = environment(.prediction_closure))
  fun
}

.as_matrix <- function(x, name = "x") {
  if (is.null(dim(x))) x <- matrix(x, ncol = 1L)
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  x
}

.check_lengths <- function(y, d, x, z) {
  n <- length(y)
  if (length(d) != n || nrow(x) != n || nrow(z) != n) {
    stop("response, diagonal, target, and instrument must have the same number of observations.", call. = FALSE)
  }
  n
}

.merge_control <- function(defaults, control) {
  if (is.null(control)) return(defaults)
  if (!is.list(control)) stop("control must be a list.", call. = FALSE)
  bad <- setdiff(names(control), names(defaults))
  if (length(bad)) stop("unknown control argument(s): ", paste(bad, collapse = ", "), call. = FALSE)
  utils::modifyList(defaults, control)
}

.safe_solve <- function(A, b, ridge = 1e-10) {
  A <- as.matrix(A)
  if (nrow(A) != ncol(A)) stop("A must be square.", call. = FALSE)
  A <- (A + t(A)) / 2
  out <- tryCatch(solve(A, b), error = function(e) NULL)
  if (!is.null(out) && all(is.finite(out))) return(out)
  solve(A + diag(ridge, nrow(A)), b)
}

.safe_inverse <- function(A, ridge = 1e-10) {
  .safe_solve(A, diag(nrow(A)), ridge = ridge)
}

.scale_fit <- function(x) {
  x <- .as_matrix(x)
  center <- colMeans(x)
  scale <- apply(x, 2L, stats::sd)
  scale[!is.finite(scale) | scale < 1e-8] <- 1
  list(center = center, scale = scale)
}

.scale_apply <- function(x, spec) {
  x <- .as_matrix(x)
  sweep(sweep(x, 2L, spec$center, "-"), 2L, spec$scale, "/")
}

.select_centers <- function(x, n_centers, seed = 1L) {
  x <- .as_matrix(x)
  k <- min(as.integer(n_centers), nrow(x))
  if (k < 1L) stop("number of centers must be positive.", call. = FALSE)
  set.seed(seed)
  x[sample.int(nrow(x), k), , drop = FALSE]
}

.median_distance <- function(x, max_n = 1500L, seed = 1L) {
  x <- .as_matrix(x)
  if (nrow(x) > max_n) {
    set.seed(seed)
    x <- x[sample.int(nrow(x), max_n), , drop = FALSE]
  }
  d <- as.numeric(stats::dist(x))
  d <- d[is.finite(d) & d > 0]
  if (!length(d)) return(1)
  stats::median(d)
}
