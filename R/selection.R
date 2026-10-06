#' Cell conditional-moment U-statistic
#'
#' Excludes self-products within each joint conditioning cell. Singleton cells
#' contribute zero; they cannot supply a pair of distinct observations.
#' @param residual Numeric residual vector or matrix, one candidate per column.
#' @param instrument Numeric conditioning matrix. All columns are retained.
#' @return Symmetric matrix, which need not be positive semidefinite.
#' @export
cell_moment_gram <- function(residual, instrument) {
  residual <- as.matrix(residual)
  instrument <- .as_matrix(instrument)
  n <- nrow(residual)
  if (!n || nrow(instrument) != n || any(!is.finite(residual)) ||
      any(!is.finite(instrument))) stop("Cell scoring requires finite, matching rows.")
  keys <- .cell_keys(instrument)
  index <- match(keys, unique(keys))
  counts <- tabulate(index)
  eligible <- counts > 1L
  if (!any(eligible)) stop("Cell U-statistic requires a cell with at least two validation observations.")
  sums <- rowsum(residual, index, reorder = FALSE)
  denominator <- as.double(n) * (counts - 1)
  # sum_{i != j} r_i r_j' = S S' - sum_i r_i r_i'.
  G <- crossprod(sums[eligible, , drop = FALSE] / sqrt(denominator[eligible]))
  rows <- eligible[index]
  G <- G - crossprod(residual[rows, , drop = FALSE] / sqrt(denominator[index[rows]]))
  G <- (G + t(G)) / 2
  attr(G, "cell_scoring") <- list(n = n, cells = length(counts),
    singleton_cells = sum(!eligible), paired_rows = sum(counts[eligible]))
  G
}

.project_gram_psd <- function(G) {
  G <- (G + t(G)) / 2
  e <- eigen(G, symmetric = TRUE)
  value <- pmax(e$values, 0)
  projected <- tcrossprod(sweep(e$vectors, 2L, sqrt(value), `*`))
  dimnames(projected) <- dimnames(G)
  list(gram = projected, eigenvalues = e$values,
       removed_negative_eigenvalues = sum(e$values < 0),
       adjustment_norm = sqrt(sum((projected - G)^2)))
}

#' Select a positive penalty on an extendable grid
#'
#' Boundary minima are recorded as failed attempts and the grid is extended.
#' Only an interior minimum is accepted. An unresolved boundary raises a
#' structured error carrying all evaluated losses and boundary attempts.
#' @param scales Positive initial penalty scales.
#' @param evaluate Function returning losses for newly requested scales.
#' @param max_extensions Maximum number of extensions before failing.
#' @param extension_factor Multiplicative width of each extension.
#' @param extension_points Number of logarithmically spaced new points per side.
#' @return List with selected scale, penalty_cv, and boundary_history.
#' @export
select_penalty_grid <- function(scales, evaluate, max_extensions = 16L,
                                extension_factor = 10, extension_points = 4L) {
  if (!length(scales) || any(!is.finite(scales)) || any(scales <= 0))
    stop("Penalty cross-validation requires positive candidate scales.")
  if (!is.function(evaluate)) stop("evaluate must be a function.")
  max_extensions <- .ensemble_integer(max_extensions, "max_extensions", 0L)
  extension_points <- .ensemble_integer(extension_points, "extension_points")
  if (length(extension_factor) != 1L || !is.finite(extension_factor) || extension_factor <= 1)
    stop("extension_factor must exceed one.")
  grid <- sort(unique(scales)); losses <- numeric(); rounds <- integer()
  events <- list(); pending <- grid
  for (attempt in 0:max_extensions) {
    current <- evaluate(pending)
    if (length(current) != length(pending) || any(!is.finite(current)))
      stop("Penalty selection requires a finite loss for every scale.")
    losses <- c(losses, current); rounds <- c(rounds, rep(attempt, length(pending)))
    evaluated <- if (attempt == 0L) pending else c(evaluated, pending)
    order <- order(evaluated); evaluated <- evaluated[order]
    losses <- losses[order]; rounds <- rounds[order]
    minima <- which(losses == min(losses))
    interior <- minima[minima > 1L & minima < length(evaluated)]
    cv <- data.frame(scale = evaluated, loss = losses, grid_round = rounds,
                     selected = FALSE)
    if (length(interior)) {
      chosen <- interior[1L]; cv$selected[chosen] <- TRUE
      return(list(scale = evaluated[chosen], penalty_cv = cv,
        boundary_history = if (length(events)) do.call(rbind, events) else data.frame(),
        status = if (attempt) "interior_after_extension" else "interior"))
    }
    low <- 1L %in% minima; high <- length(evaluated) %in% minima
    events[[length(events) + 1L]] <- data.frame(attempt = attempt,
      lower = min(evaluated), upper = max(evaluated), minimum_loss = min(losses),
      lower_boundary = low, upper_boundary = high, status = "boundary_failure")
    if (attempt == max_extensions) {
      stop(structure(list(message = paste0("Penalty grid has a boundary minimum after ",
        attempt, " extensions; no interior penalty was accepted."), call = NULL,
        penalty_cv = cv, boundary_history = do.call(rbind, events)),
        class = c("cmbridge_penalty_boundary_error", "error", "condition")))
    }
    step <- seq_len(extension_points) / extension_points
    pending <- c(if (low) min(evaluated) / extension_factor^step,
                 if (high) max(evaluated) * extension_factor^step)
    pending <- sort(setdiff(pending, evaluated))
    if (!length(pending) || any(!is.finite(pending)) || any(pending <= 0))
      stop("Penalty grid extension exceeded the positive finite numerical range.")
  }
}
