# Parameter-choice rules for the empirical inverse problem. No population
# functions or validation rows enter this path; ensemble stacking stays RKHS.

.inverse_selection_control <- function(control) {
  out <- .merge_control(list(scales = 10^seq(-10, 0, by = .25),
    iteration_grid = unique(round(exp(seq(log(2), log(4000), length.out = 25L)))),
    max_extensions = 6L, extension_factor = 10, extension_points = 4L,
    max_iterations = 50000L, derivative_window = 5L,
    trace_probes = 8L, trace_step = 1e-4, seed = 1L), control)
  out$max_extensions <- .ensemble_integer(out$max_extensions, "max_extensions", 0L)
  out$extension_points <- .ensemble_integer(out$extension_points, "extension_points", 1L)
  out$derivative_window <- .ensemble_integer(out$derivative_window, "derivative_window", 3L)
  out$trace_probes <- .ensemble_integer(out$trace_probes, "trace_probes", 1L)
  out$max_iterations <- .ensemble_integer(out$max_iterations, "max_iterations", 2L)
  if (!is.finite(out$extension_factor) || out$extension_factor <= 1 ||
      !is.finite(out$trace_step) || out$trace_step <= 0)
    stop("Invalid inverse selection extension or trace controls.")
  out
}

.inverse_native_control <- function(control, method, scale) {
  for (name in c("penalty_selection", "penalty_selection_control", "penalty_scales",
    "penalty_ids", "penalty_folds", "penalty_max_extensions", "penalty_extension_factor",
    "penalty_extension_points", "iteration_grid")) control[[name]] <- NULL
  if (method == "landweber") {
    control$n_iter <- max(1L, as.integer(round(1 / scale)))
    control$max_iter <- max(2000L, control$n_iter)
  } else control$lambda <- scale
  control
}

# Orthogonal compression removes redundant critic coordinates. It preserves
# the native weighted moment norm, rather than dividing by redundant columns.
.inverse_geometry <- function(model, y, d, x, z) {
  n <- length(y); active <- d != 0
  Q <- model$instrument_features(z); W <- model$moment_weight
  ew <- eigen((W + t(W)) / 2, symmetric = TRUE)
  use <- ew$values > max(1e-14, max(ew$values) * 1e-12)
  root <- sweep(t(ew$vectors[, use, drop = FALSE]), 1L, sqrt(ew$values[use]), `*`)
  L <- root %*% t(Q) / n
  el <- eigen(tcrossprod(L), symmetric = TRUE)
  use <- el$values > max(1e-20, max(el$values) * 1e-11)
  if (!any(use)) stop("The empirical critic has zero rank.")
  L <- crossprod(el$vectors[, use, drop = FALSE], L)
  cell <- model$method == "sieve_md" && identical(model$target_spec$type, "cell") &&
    (model$tuning$link == "identity" || model$tuning$penalty_on == "function") &&
    model$target_design_cols == length(model$target_spec$levels) + 1L
  if (model$method == "pmmr") {
    H <- .nystrom_apply(.scale_apply(x[active, , drop = FALSE], model$target_scale), model$target_nystrom)
  } else if (cell) {
    H <- .eval_basis_spec(x[active, , drop = FALSE], model$target_spec)[, -1L, drop = FALSE]
  } else H <- .eval_basis_spec(x[active, , drop = FALSE], model$target_spec)
  theta <- if (cell) {
    if (model$tuning$link == "identity") model$coefficients[1L] + model$coefficients[-1L] else model$solver$parameters
  } else model$coefficients
  eta <- as.numeric(H %*% theta)
  link <- model$tuning$link
  values <- if (link == "identity") eta else .cm_link(eta, link)
  first <- if (link == "identity") rep(1, length(eta)) else .cm_link_derivative(eta, link)
  second <- if (link == "identity") rep(0, length(eta)) else exp(if (link == "inverse_logit") -eta else eta)
  P <- diag(length(theta))
  if (cell) P <- P - matrix(1 / length(theta), length(theta), length(theta)) else
    if (model$method == "sieve_md") {
      if (link != "identity") {
        scales <- sqrt(colMeans(H^2)); scales[scales < 1e-8] <- 1; scales[1L] <- 1
        weights <- rep(1, ncol(H))
        if (model$target_spec$type == "cell_linear") {
          categories <- seq.int(model$target_spec$p + 2L, ncol(H))
          scales[categories] <- 1; weights[categories] <- model$tuning$cell_penalty_multiplier
        }
        P <- diag(weights * scales^2)
      }
      P[1L, 1L] <- 0
    }
  prediction <- numeric(n); prediction[active] <- values
  residual <- y - d * prediction
  e <- as.numeric(L %*% residual)
  J <- L[, active, drop = FALSE] %*% (d[active] * first * H)
  S <- crossprod(J)
  if (cell && link != "identity") {
    cell_values <- .cm_link(theta, link)
    cell_first <- .cm_link_derivative(theta, link)
    cell_second <- exp(if (link == "inverse_logit") -theta else theta)
    penalty_values <- as.numeric(P %*% cell_values)
    derivative <- diag(cell_first, length(cell_first))
    penalty_hessian <- derivative %*% P %*% derivative +
      diag(cell_second * penalty_values, length(cell_second))
    seminorm2 <- as.numeric(crossprod(cell_values, P %*% cell_values))
  } else {
    penalty_hessian <- P
    seminorm2 <- as.numeric(crossprod(theta, P %*% theta))
  }
  if (model$method == "landweber") {
    # Landweber regularizes deviations from its initialized coefficient vector.
    transform <- diag(ncol(H))
    if (link != "identity" && isTRUE(model$tuning$precondition)) {
      ee <- eigen(crossprod(H) / nrow(H), symmetric = TRUE)
      transform <- ee$vectors %*% diag(1 / sqrt(pmax(ee$values, max(ee$values) * 1e-3))) %*% t(ee$vectors)
    }
    initial <- numeric(ncol(H))
    if (isTRUE(model$tuning$init_intercept)) {
      value <- mean(y) / mean(d)
      initial[1L] <- if (link == "identity") value else {
        gap <- max(.1, value - if (link == "inverse_logit") 1 else 0)
        if (link == "inverse_logit") -log(gap) else log(gap)
      }
    }
    seminorm2 <- sum(solve(transform, theta - initial)^2)
  }
  list(L = L, H = H, J = J, S = S, first = first, second = second,
    theta = theta, P = P, penalty_hessian = penalty_hessian, active = active,
    residual = residual, e = e, values = values, seminorm = sqrt(max(0, seminorm2)),
    rank = nrow(L), native_loss = sum(e^2), link = link)
}

.inverse_ridge_df <- function(g, model, d) {
  curvature <- g$S - crossprod(g$H,
    d[g$active] * g$second * as.numeric(crossprod(g$L[, g$active, drop = FALSE], g$e)) * g$H) +
      model$tuning$lambda * g$penalty_hessian
  ee <- eigen((curvature + t(curvature)) / 2, symmetric = TRUE)
  tolerance <- max(1e-14, max(abs(ee$values)) * 1e-10)
  if (min(ee$values) < -tolerance) stop("The fitted local Hessian is not positive semidefinite.")
  use <- ee$values > tolerance
  sum(colSums((g$J %*% ee$vectors[, use, drop = FALSE])^2) / ee$values[use])
}

.inverse_landweber_df <- function(g, model, y, d, x, z, control, scale, options) {
  if (g$link == "identity") {
    ee <- eigen((g$S + t(g$S)) / 2, symmetric = TRUE)
    s <- pmax(ee$values, 0); step <- model$tuning$step
    remainder <- (1 - step * s)^model$tuning$iterations_used
    df <- sum(1 - remainder)
    if (isTRUE(model$tuning$init_intercept)) {
      # The initialized intercept depends on the moment response as well.
      v <- .safe_solve(tcrossprod(g$L), as.numeric(g$L %*% rep(1 / length(y), length(y))))
      initial <- numeric(ncol(g$H)); initial[1L] <- 1 / mean(d)
      propagated <- ee$vectors %*% (remainder * as.numeric(crossprod(ee$vectors, initial)))
      df <- df + sum(v * as.numeric(g$J %*% propagated))
    }
    return(list(df = df, kind = "exact_fixed_iteration", trace_se = 0))
  }
  # Differentiate the actual nonlinear iteration, including initialization,
  # backtracking and preconditioning. Common probes are used along the path.
  set.seed(options$seed)
  probes <- if (g$rank <= options$trace_probes) diag(g$rank) else
    matrix(sample(c(-1, 1), g$rank * options$trace_probes, replace = TRUE), g$rank)
  lift <- crossprod(g$L, .safe_solve(tcrossprod(g$L), probes))
  epsilon <- options$trace_step * max(1, sqrt(sum((g$L %*% y)^2) / g$rank))
  estimates <- vapply(seq_len(ncol(probes)), function(j) {
    yy <- lapply(c(-1, 1), function(sign) y + sign * epsilon * lift[, j])
    fits <- lapply(yy, function(response) fit_cm(response, d, x, z, "landweber",
      .inverse_native_control(control, "landweber", scale)))
    delta <- as.numeric(g$L %*% (d * (ifelse(is.na(fits[[2L]]$fitted), 0, fits[[2L]]$fitted) -
      ifelse(is.na(fits[[1L]]$fitted), 0, fits[[1L]]$fitted)))) / (2 * epsilon)
    sum(probes[, j] * delta)
  }, 0)
  exact <- g$rank <= options$trace_probes
  list(df = if (exact) sum(estimates) else mean(estimates),
    kind = if (exact) "central_difference_trace" else "common_probe_trace",
    trace_se = if (exact || length(estimates) < 2L) 0 else stats::sd(estimates) / sqrt(length(estimates)))
}

.inverse_derivatives <- function(parameter, values, width) {
  t <- log(parameter); values <- as.matrix(values)
  first <- second <- matrix(NA_real_, nrow(values), ncol(values))
  for (j in seq_along(t)) {
    indices <- order(abs(t - t[j]))[seq_len(min(width, length(t)))]
    offset <- t[indices] - t[j]
    design <- cbind(1, offset, offset^2)
    coefficient <- qr.solve(design, values[indices, , drop = FALSE])
    first[j, ] <- coefficient[2L, ]; second[j, ] <- 2 * coefficient[3L, ]
  }
  list(first = first, second = second)
}

.inverse_path_scores <- function(grid, records, rule, options) {
  valid <- vapply(records, function(q) is.finite(q$residual_norm) && is.finite(q$seminorm), TRUE)
  score <- rep(Inf, length(grid)); curvature <- change <- rep(NA_real_, length(grid))
  if (rule == "gcv") {
    score <- vapply(records, `[[`, 0, "gcv")
  } else if (sum(valid) >= options$derivative_window) {
    at <- which(valid)
    if (rule == "lcurve") {
      at <- at[vapply(records[at], function(q) q$residual_norm > 1e-14 && q$seminorm > 1e-14, TRUE)]
      if (length(at) >= options$derivative_window) {
        points <- cbind(log(vapply(records[at], `[[`, 0, "residual_norm")),
          log(vapply(records[at], `[[`, 0, "seminorm")))
        derivative <- .inverse_derivatives(grid[at], points, options$derivative_window)
        dx <- derivative$first[, 1L]; dy <- derivative$first[, 2L]
        ddx <- derivative$second[, 1L]; ddy <- derivative$second[, 2L]
        # Increasing lambda traverses this curve clockwise; positive corner
        # curvature therefore has the opposite sign of the usual determinant.
        curvature[at] <- (dy * ddx - dx * ddy) / pmax((dx^2 + dy^2)^1.5, 1e-30)
        score[at] <- ifelse(is.finite(curvature[at]) & curvature[at] > 0, -curvature[at], Inf)
      }
    } else {
      functions <- do.call(rbind, lapply(records[at], `[[`, "values"))
      derivative <- .inverse_derivatives(grid[at], functions, options$derivative_window)$first
      change[at] <- sqrt(rowMeans(derivative^2))
      # The criterion is the RMS derivative of the fitted function in log
      # regularization. Only strict interior valleys are admissible: neither
      # a constant oversmoothed limit nor numerical convergence counts as one.
      for (j in seq.int(2L, length(at) - 1L)) {
        a <- at[j]; value <- change[a]
        noise_floor <- 1e-8 * max(1, sqrt(mean(functions[j, ]^2)))
        if (value > noise_floor && value < min(change[at[j - 1L]], change[at[j + 1L]]) * (1 - 1e-3)) score[a] <- value
      }
    }
  }
  list(score = score, curvature = curvature, change = change)
}

#' Select regularization for a conditional-moment learner
#'
#' Compare generalized cross-validation, L-curve curvature, or interior
#' quasi-optimality along a positive training-only regularization path.
#' @param response Numeric conditional-moment loading.
#' @param diagonal Numeric multiplier, or NULL for ones.
#' @param target,instrument Target and conditioning predictor matrices.
#' @param method One of sieve_md, landweber, pmmr.
#' @param control Native learner controls. penalty_scales supplies a ridge grid.
#' @param rule One of gcv, lcurve, quasiopt.
#' @param selection_control Path, derivative and nonlinear trace controls.
#' @return A cmbridge_fit with the complete penalty_path, selection diagnostics
#'   and recorded boundary failures. An unresolved boundary raises an error.
#' @export
select_cm_penalty <- function(response, diagonal = NULL, target, instrument,
                             method = c("sieve_md", "landweber", "pmmr"), control = list(),
                             rule = c("gcv", "lcurve", "quasiopt"), selection_control = list()) {
  method <- match.arg(method); rule <- match.arg(rule)
  y <- as.numeric(response); d <- if (is.null(diagonal)) rep(1, length(y)) else as.numeric(diagonal)
  x <- .as_matrix(target); z <- .as_matrix(instrument); .check_lengths(y, d, x, z)
  options <- .inverse_selection_control(selection_control)
  saved <- if (exists(".Random.seed", .GlobalEnv)) get(".Random.seed", .GlobalEnv) else NULL
  on.exit(if (is.null(saved)) {if (exists(".Random.seed", .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)}
    else assign(".Random.seed", saved, .GlobalEnv), add = TRUE)
  if (method == "landweber") {
    iterations <- if (is.null(control$iteration_grid)) options$iteration_grid else control$iteration_grid
    if (any(!is.finite(iterations)) || any(iterations < 1) || any(iterations != round(iterations))) stop("Invalid iteration grid.")
    grid <- sort(unique(1 / iterations))
  } else grid <- if (is.null(control$penalty_scales)) options$scales else control$penalty_scales
  if (length(unique(grid)) < 5L || any(!is.finite(grid)) || any(grid <= 0)) stop("Use at least five distinct positive path scales.")
  grid <- sort(unique(grid)); cache <- new.env(parent = emptyenv()); boundaries <- list()
  fit_at <- function(scale) {
    key <- format(scale, digits = 17)
    if (exists(key, cache, inherits = FALSE)) return(get(key, cache))
    result <- tryCatch({
      model <- fit_cm(y, d, x, z, method, .inverse_native_control(control, method, scale))
      g <- .inverse_geometry(model, y, d, x, z)
      stopifnot(abs(g$native_loss - model$moment_loss) < 1e-8 * max(1, model$moment_loss))
      trace <- if (rule != "gcv") list(df = NA_real_, kind = "not_required", trace_se = NA_real_) else
        tryCatch(if (method == "landweber") .inverse_landweber_df(g, model, y, d, x, z, control, scale, options) else
          list(df = .inverse_ridge_df(g, model, d), kind = if (g$link == "identity") "exact_ridge" else "local_hessian", trace_se = 0),
          error = function(e) list(df = NA_real_, kind = conditionMessage(e), trace_se = NA_real_))
      denominator <- 1 - trace$df / g$rank
      gcv <- if (is.finite(denominator) && denominator > 1e-8 && trace$df >= -1e-8)
        g$native_loss / g$rank / denominator^2 else Inf
      list(model = model, residual_norm = sqrt(g$native_loss), seminorm = g$seminorm,
        values = g$values, rank = g$rank, df = trace$df, trace_kind = trace$kind,
        trace_se = trace$trace_se, gcv = gcv, failure = "")
    }, error = function(e) list(model = NULL, residual_norm = Inf, seminorm = Inf,
      values = numeric(), rank = NA_integer_, df = NA_real_, trace_kind = "failed",
      trace_se = NA_real_, gcv = Inf, failure = conditionMessage(e)))
    assign(key, result, cache); result
  }
  for (attempt in 0:options$max_extensions) {
    records <- lapply(grid, fit_at); scores <- .inverse_path_scores(grid, records, rule, options)
    path <- data.frame(scale = grid, score = scores$score,
      residual_norm = vapply(records, `[[`, 0, "residual_norm"),
      seminorm = vapply(records, `[[`, 0, "seminorm"),
      degrees_of_freedom = vapply(records, `[[`, 0, "df"),
      critic_rank = vapply(records, `[[`, 0L, "rank"),
      trace_kind = vapply(records, `[[`, "", "trace_kind"),
      trace_se = vapply(records, `[[`, 0, "trace_se"),
      curvature = scores$curvature, change = scores$change,
      failed = vapply(records, function(q) is.null(q$model), TRUE),
      failure = vapply(records, `[[`, "", "failure"), selected = FALSE)
    choice <- if (any(is.finite(scores$score))) which.min(scores$score) else NA_integer_
    accepted <- !is.na(choice) && choice > 1L && choice < length(grid) &&
      !any(path$failed[c(choice - 1L, choice + 1L)])
    if (accepted) {
      path$selected[choice] <- TRUE; model <- records[[choice]]$model
      model$penalty_selection <- rule; model$penalty_path <- path
      model$penalty_cv <- path; model$penalty_cv$loss <- path$score
      model$penalty_boundary <- if (length(boundaries)) do.call(rbind, boundaries) else data.frame()
      model$penalty_status <- "interior_inverse_parameter_choice"
      model$penalty_parameter <- if (method == "landweber") "inverse_iterations" else "lambda"
      model$penalty_selection_control <- options; model$penalty_training_n <- length(y)
      model$call <- match.call(); return(model)
    }
    boundaries[[length(boundaries) + 1L]] <- data.frame(attempt = attempt, lower = min(grid), upper = max(grid),
      selected_scale = if (is.na(choice)) NA_real_ else grid[choice],
      status = if (is.na(choice)) "no_interior_corner_or_valley" else "boundary_or_failed_neighbor")
    if (attempt == options$max_extensions) break
    step <- seq_len(options$extension_points) / options$extension_points
    low <- min(grid) / options$extension_factor^step
    high <- max(grid) * options$extension_factor^step
    added <- if (is.na(choice)) c(low, high) else if (choice == 1L) low else if (choice == length(grid)) high else c(low, high)
    if (method == "landweber") added <- 1 / pmin(options$max_iterations, pmax(1, round(1 / added)))
    next_grid <- sort(unique(c(grid, added)))
    if (identical(next_grid, grid)) break
    grid <- next_grid
  }
  cv <- path; cv$loss <- path$score
  stop(structure(list(message = paste("No admissible interior", rule, "parameter choice."), call = NULL,
    penalty_path = path, penalty_cv = cv,
    boundary_history = do.call(rbind, boundaries)), class = c("cmbridge_penalty_selection_error", "error", "condition")))
}
