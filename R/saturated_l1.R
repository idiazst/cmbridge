# A categorical conditional-moment model. Saturation refers to joint cells,
# rather than marginal polynomial terms. The penalty is on cell deviations
# from an unpenalized constant function.
.cell_keys <- function(x, columns = seq_len(ncol(x))) {
  if (!length(columns)) return(rep("constant", nrow(x)))
  x <- x[, columns, drop = FALSE]
  # Convert each distinct discrete value once. Numeric-to-text conversion
  # on every repeated row dominated large finite-support study runtimes.
  columns <- lapply(as.data.frame(x), function(column) {
    values <- unique(column)
    as.character(values)[match(column, values)]
  })
  do.call(paste, c(columns, sep = "|"))
}

.fit_saturated_l1 <- function(y, d, x, z, control) {
  ctrl <- .merge_control(list(penalty = 0.001, target_columns = seq_len(ncol(x)),
                             instrument_columns = seq_len(ncol(z)),
                             lower = -Inf, upper = Inf, tolerance = 1e-8,
                             max_iter = 1000000L), control)
  if (length(ctrl$penalty) != 1L || !is.finite(ctrl$penalty) || ctrl$penalty < 0) {
    stop("penalty must be a nonnegative finite scalar.", call. = FALSE)
  }
  if (length(ctrl$lower) != 1L || length(ctrl$upper) != 1L ||
      is.na(ctrl$lower) || is.na(ctrl$upper) || ctrl$lower > ctrl$upper) {
    stop("lower and upper must be ordered scalar prediction bounds.", call. = FALSE)
  }
  active <- which(d != 0)
  if (!length(active)) stop("At least one nonzero diagonal is required.")
  xkey <- .cell_keys(x[active, , drop = FALSE], ctrl$target_columns)
  zkey <- .cell_keys(z, ctrl$instrument_columns)
  xlevels <- sort(unique(xkey))
  zlevels <- sort(unique(zkey))
  zi <- match(zkey, zlevels)
  xi <- match(xkey, xlevels)
  counts <- tabulate(zi, length(zlevels))
  sums <- Matrix::sparseMatrix(i = zi[active], j = xi, x = d[active],
                              dims = c(length(zlevels), length(xlevels)))
  A <- Matrix::Diagonal(x = 1 / counts) %*% sums
  cvec <- as.numeric(rowsum(y, zi, reorder = FALSE)[match(seq_along(zlevels), unique(zi)), 1]) / counts
  design <- cbind(Matrix::Matrix(as.numeric(Matrix::rowSums(A)), ncol = 1, sparse = TRUE), A)
  # A constant h is multiplied by D in the moment equation, so glmnet's
  # ordinary intercept would be incorrect when D is an observation indicator.
  if (all(cvec == 0)) {
    coefficient <- rep(0, length(xlevels) + 1L)
  } else if (nrow(design) == 1L) {
    coefficient <- c(if (design[1, 1] != 0) cvec[1] / design[1, 1] else 0,
                     rep(0, length(xlevels)))
  } else {
    # A decreasing path supplies warm starts for small penalties. Fitting one
    # tiny lambda from zero is particularly slow for nearly singular operators.
    path <- sort(unique(c(0.1, 0.01, 0.001, ctrl$penalty)), decreasing = TRUE)
    path <- path[path >= ctrl$penalty]
    fit <- glmnet::glmnet(design, cvec, weights = counts / sum(counts),
                         family = "gaussian", intercept = FALSE, standardize = FALSE,
                         penalty.factor = c(0, rep(1, length(xlevels))),
                         lambda = path, thresh = ctrl$tolerance,
                         maxit = ctrl$max_iter)
    if (fit$jerr != 0 || !length(fit$lambda) ||
        abs(tail(fit$lambda, 1) - ctrl$penalty) > 1e-10) {
      stop("Saturated L1 conditional-moment solver failed to converge at the requested penalty.", call. = FALSE)
    }
    coefficient <- as.numeric(fit$beta[, length(fit$lambda)])
  }
  predict_fun <- .prediction_closure(function(newx) {
    newx <- .as_matrix(newx)
    index <- match(.cell_keys(newx, ctrl$target_columns), xlevels)
    deviation <- numeric(nrow(newx))
    seen <- !is.na(index)
    deviation[seen] <- coefficient[1L + index[seen]]
    coefficient[1L] + deviation
  }, list(ctrl = ctrl, xlevels = xlevels, coefficient = coefficient))
  pred <- rep(NA_real_, length(y))
  pred[active] <- predict_fun(x[active, , drop = FALSE])
  residual <- y - d * ifelse(is.na(pred), 0, pred)
  means <- as.numeric(rowsum(residual, zi, reorder = FALSE)[match(seq_along(zlevels), unique(zi)), 1]) / counts
  list(coefficients = coefficient, target_levels = xlevels, instrument_levels = zlevels,
       target_spec = list(type = "cells"), fitted = pred, residual = residual,
       target_design_cols = length(xlevels), moment_loss = sum(counts * means^2) / length(y),
       tuning = ctrl, predict_fun = predict_fun,
       instrument_features = .prediction_closure(function(newz) {
         key <- match(.cell_keys(.as_matrix(newz), ctrl$instrument_columns), zlevels)
         keep <- which(!is.na(key))
         Matrix::sparseMatrix(i = keep, j = key[keep], x = 1,
                             dims = c(nrow(newz), length(zlevels)))
       }, list(ctrl = ctrl, zlevels = zlevels)))
}
