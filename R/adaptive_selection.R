# Experimental adaptations of Ghassami--Robins--Rotnitzky (2025, Eq. 6)
# and Tan--Syrgkanis (2026, discrepancy principle). Native learner objectives
# are retained. No population truth is consulted by either selector.

.adaptive_control <- function(control) {
  out <- .merge_control(list(projection_basis = "cell_linear", projection_degree = 1L,
    projection_df = 6L, projection_ridges = 10^seq(-6, 1, length.out = 15L),
    critic_ridge = 1e-6, bootstrap_draws = 199L, noise_probability = .95,
    max_extensions = 16L, extension_factor = 10, extension_points = 4L,
    iteration_grid = unique(round(exp(seq(log(2), log(2000), length.out = 21L)))),
    max_iterations = 100000L, limit_tolerance = 1e-8, seed = 1L), control)
  out$projection_basis <- match.arg(out$projection_basis, c("cell_linear", "cell", "poly", "bs"))
  if (any(!is.finite(out$projection_ridges)) || any(out$projection_ridges <= 0))
    stop("Projection regression requires positive ridge values.")
  if (!is.finite(out$critic_ridge) || out$critic_ridge <= 0 ||
      !is.finite(out$noise_probability) || out$noise_probability <= 0 || out$noise_probability >= 1)
    stop("Invalid critic ridge or noise_probability.")
  out$bootstrap_draws <- .ensemble_integer(out$bootstrap_draws, "bootstrap_draws", 19L)
  out$seed <- .ensemble_integer(out$seed, "selection seed", 0L)
  out
}

.adaptive_labels <- function(ids, groups = 3L) {
  values <- as.character(ids)
  codes <- suppressWarnings(abs(as.numeric(values)))
  for (j in which(!is.finite(codes))) codes[j] <- strtoi(substr(
    digest::digest(values[j], algo = "xxhash32"), 1L, 7L), 16L)
  # Stable person labels agree on shared rows, including after restriction to
  # measured observations. Use the first digit varying in a deeper sample.
  for (digit in 0:32) {
    labels <- as.integer(floor(codes / groups^digit) %% groups) + 1L
    if (length(unique(labels)) >= min(groups, length(unique(values)))) return(labels)
  }
  stop("Could not construct distinct person folds.")
}

.adaptive_plans <- function(index, ids, labels = NULL) {
  if (is.null(labels) || length(unique(labels)) < 2L) labels <- .adaptive_labels(ids[index])
  if (length(labels) != length(index) || anyNA(labels)) stop("Invalid adaptive selection folds.")
  if (any(labels != labels[match(as.character(ids[index]), as.character(ids[index]))]))
    stop("A person occurs in more than one adaptive selection fold.")
  levels <- sort(unique(labels)); plans <- list()
  for (j in seq_along(levels)) {
    validation <- index[labels == levels[j]]
    if (length(levels) >= 3L) {
      projection <- index[labels == levels[j %% length(levels) + 1L]]
      train <- setdiff(index, c(validation, projection))
      plans[[length(plans) + 1L]] <- list(train = train, projection = projection, validation = validation)
    } else {
      train <- index[labels != levels[j]]
      split <- .adaptive_labels(ids[validation], 2L)
      for (side in 1:2) plans[[length(plans) + 1L]] <- list(train = train,
        projection = validation[split == side], validation = validation[split != side])
    }
  }
  for (plan in plans) {
    if (any(lengths(plan) < 2L)) stop("Each adaptive selection sample needs at least two rows.")
    people <- lapply(plan, function(rows) unique(as.character(ids[rows])))
    if (length(intersect(people$train, people$projection)) ||
        length(intersect(people$train, people$validation)) ||
        length(intersect(people$projection, people$validation))) stop("Adaptive selection samples overlap.")
  }
  plans
}

.adaptive_design <- function(z, newz, control) {
  spec <- .fit_basis_spec(z, control$projection_basis, control$projection_degree, control$projection_df)
  train <- .eval_basis_spec(z, spec); validation <- .eval_basis_spec(newz, spec)
  Q <- train[, -1L, drop = FALSE]; N <- validation[, -1L, drop = FALSE]
  center <- colMeans(Q); Q <- sweep(Q, 2L, center); N <- sweep(N, 2L, center)
  scales <- sqrt(colMeans(Q^2)); scales[scales < 1e-12] <- 1
  Q <- sweep(Q, 2L, scales, `/`); N <- sweep(N, 2L, scales, `/`)
  e <- eigen(crossprod(Q) / nrow(Q), symmetric = TRUE)
  keep <- e$values > max(1e-12, max(e$values) * 1e-12)
  U <- Q %*% e$vectors[, keep, drop = FALSE]; NU <- N %*% e$vectors[, keep, drop = FALSE]
  critic <- eigen(crossprod(train) / nrow(train), symmetric = TRUE)
  retained <- critic$values > max(1e-12, max(critic$values) * 1e-12)
  transform <- sweep(critic$vectors[, retained, drop = FALSE], 2L,
    sqrt(critic$values[retained] + control$critic_ridge), `/`)
  list(train = train, validation = validation, spec = spec,
    U = U, NU = NU, eigenvalues = pmax(e$values[keep], 0),
    critic_features = validation %*% transform, critic_rank = sum(retained))
}

# Ordinary leave-one-out ridge prediction loss on the projection-training
# sample. The scoring sample is never used to learn features or select ridge.
.adaptive_projection <- function(design, response, ridges) {
  response <- as.matrix(response); n <- nrow(response)
  mean_response <- colMeans(response); Y <- sweep(response, 2L, mean_response)
  variance <- colMeans(Y^2); variance[variance < 1e-12] <- 1
  if (!ncol(design$U)) return(list(predictions = matrix(mean_response, nrow(design$validation), ncol(Y), byrow = TRUE),
    ridge = ridges[1L], cv_loss = 0))
  U <- design$U; NU <- design$NU
  rhs <- crossprod(U, Y) / n
  cache <- new.env(parent = emptyenv())
  evaluate <- function(requested) vapply(requested, function(ridge) {
    denominator <- design$eigenvalues + ridge
    coefficients <- sweep(rhs, 1L, denominator, `/`)
    predicted <- U %*% coefficients
    leverage <- 1 / n + rowSums(sweep(U^2, 2L, denominator, `/`)) / n
    residual <- sweep(Y - predicted, 1L, pmax(1 - leverage, 1e-10), `/`)
    value <- list(loss = mean(sweep(residual^2, 2L, variance, `/`)), coefficients = coefficients)
    assign(format(ridge, digits = 17), value, cache)
    value$loss
  }, 0)
  selected <- select_penalty_grid(ridges, evaluate)
  value <- get(format(selected$scale, digits = 17), cache)
  list(predictions = sweep(NU %*% value$coefficients, 2L, mean_response, `+`),
    ridge = selected$scale, cv_loss = value$loss, ridge_boundary = FALSE,
    ridge_status = selected$status)
}

.adaptive_values <- function(model, x, d, rows) {
  answer <- numeric(length(rows)); active <- which(d[rows] != 0)
  answer[active] <- d[rows[active]] * predict(model, x[rows[active], , drop = FALSE])
  if (any(!is.finite(answer))) stop("Nonfinite adaptive selection predictions.")
  answer
}

.adaptive_score <- function(up, uv, yp, yv, zp, zv, control, known_response = FALSE,
                            design = NULL, criterion = "robins") {
  saved <- if (exists(".Random.seed", .GlobalEnv)) get(".Random.seed", .GlobalEnv) else NULL
  on.exit(if (is.null(saved)) {if (exists(".Random.seed", .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)}
    else assign(".Random.seed", saved, .GlobalEnv), add = TRUE)
  set.seed(control$seed)
  up <- as.matrix(up); uv <- as.matrix(uv)
  if (is.null(design)) design <- .adaptive_design(zp, zv, control)
  p <- ncol(up); E <- uv - yv; n <- length(yv)
  projected <- NULL; constant <- 0
  if (criterion == "robins") {
    projected <- .adaptive_projection(design, if (known_response) up else cbind(up, yp), control$projection_ridges)
    t <- projected$predictions[, seq_len(p), drop = FALSE]
    r <- if (known_response) rep(yp[1L], length(yv)) else projected$predictions[, p + 1L]
    F <- t - r; constant <- mean((r - yv)^2)
  } else F <- matrix(0, n, p)
  # Eq. 6 minus the common, candidate-independent E[(rhat - Y)^2].
  robins <- (crossprod(F, E) + crossprod(E, F) - crossprod(F)) / n
  Q <- design$critic_features
  moment <- crossprod(Q, E) / n
  discrepancy <- crossprod(moment)
  # A simultaneous multiplier-bootstrap noise envelope over all columns.
  multipliers <- matrix(stats::rnorm(n * control$bootstrap_draws), n, control$bootstrap_draws)
  noise <- matrix(0, control$bootstrap_draws, p)
  for (j in seq_len(p)) {
    contribution <- Q * E[, j]
    centered <- sweep(contribution, 2L, colMeans(contribution))
    noise[, j] <- colSums((crossprod(centered, multipliers) / n)^2)
  }
  list(robins = (robins + t(robins)) / 2, discrepancy = discrepancy,
    norm = crossprod(uv) / n, noise = noise, n = n, projection = projected[c("ridge", "cv_loss", "ridge_boundary")],
    projection_constant = constant, critic_rank = design$critic_rank)
}

.adaptive_fit_at <- function(y, d, x, z, rows, candidate, scale, control) {
  native <- candidate$control
  for (name in c("penalty_scales", "penalty_ids", "penalty_folds", "penalty_max_extensions",
    "penalty_extension_factor", "penalty_extension_points", "iteration_grid")) native[[name]] <- NULL
  if (candidate$method == "landweber") {
    iterations <- max(1L, as.integer(round(1 / scale)))
    if (!is.finite(iterations) || iterations > control$max_iterations) stop("Iteration grid exceeded max_iterations.")
    native$n_iter <- iterations; native$max_iter <- max(iterations, 2000L)
  } else if (candidate$method %in% c("sieve_md", "pmmr")) native$lambda <- scale else
    stop("Adaptive selection supports sieve_md, landweber and pmmr.")
  fit_cm(y, d[rows], x[rows, , drop = FALSE], z[rows, , drop = FALSE], candidate$method, native)
}

.adaptive_grid <- function(candidate, control) {
  if (candidate$method == "landweber") {
    iterations <- candidate$control$iteration_grid
    if (is.null(iterations)) iterations <- control$iteration_grid
    if (any(!is.finite(iterations)) || any(iterations < 1) || any(iterations != round(iterations)))
      stop("iteration_grid must contain positive integers.")
    return(sort(unique(1 / iterations)))
  }
  grid <- candidate$control$penalty_scales
  if (is.null(grid)) grid <- 10^seq(-10, 0, by = .25)
  if (any(!is.finite(grid)) || any(grid <= 0)) stop("Adaptive penalties must be positive.")
  sort(unique(grid))
}

.adaptive_discrepancy_grid <- function(grid, evaluate, cache, control, method) {
  boundaries <- list(); previous_limit <- NULL
  for (round in 0:control$max_extensions) {
    evaluate(grid[!vapply(grid, function(s) exists(format(s, digits = 17), cache, inherits = FALSE), TRUE)])
    values <- lapply(grid, function(s) get(format(s, digits = 17), cache, inherits = FALSE))
    good <- vapply(values, function(v) is.finite(v$loss), TRUE)
    if (!any(good)) stop("All discrepancy penalty trials failed.")
    noise <- do.call(cbind, lapply(values[good], `[[`, "noise"))
    threshold <- as.numeric(stats::quantile(apply(noise, 1L, max), control$noise_probability, names = FALSE))
    losses <- vapply(values, `[[`, 0, "loss")
    accepted <- which(good & losses <= threshold)
    choice <- if (length(accepted)) max(accepted) else NA_integer_
    limit <- !is.na(choice) && choice == length(grid) && (
      (method == "landweber" && grid[choice] >= 1) ||
      (!is.null(previous_limit) && abs(losses[choice] - previous_limit) <=
        control$limit_tolerance * max(1, abs(previous_limit))))
    if (!is.na(choice) && (choice < length(grid) || limit)) {
      cv <- data.frame(scale = grid, loss = losses, noise_threshold = threshold,
        admissible = good & losses <= threshold, selected = seq_along(grid) == choice,
        failed = !good, failure = vapply(values, `[[`, "", "failure"))
      return(list(scale = grid[choice], penalty_cv = cv,
        status = if (limit) "discrepancy_regularized_limit" else "discrepancy_crossing",
        boundary_history = if (length(boundaries)) do.call(rbind, boundaries) else data.frame(), threshold = threshold))
    }
    boundaries[[length(boundaries) + 1L]] <- data.frame(attempt = round, lower = min(grid),
      upper = max(grid), threshold = threshold, status = if (is.na(choice)) "no_admissible_penalty" else "unbracketed_threshold")
    if (round == control$max_extensions) stop(structure(list(message = "Discrepancy rule did not bracket its noise threshold.",
      call = NULL, boundary_history = do.call(rbind, boundaries),
      penalty_cv = data.frame(scale = grid, loss = losses, noise_threshold = threshold,
        admissible = good & losses <= threshold, selected = FALSE, failed = !good,
        failure = vapply(values, `[[`, "", "failure")), threshold = threshold),
      class = c("cmbridge_discrepancy_error", "error", "condition")))
    if (!is.na(choice)) previous_limit <- losses[choice]
    step <- seq_len(control$extension_points) / control$extension_points
    added <- if (is.na(choice)) min(grid) / control$extension_factor^step else max(grid) * control$extension_factor^step
    if (method == "landweber") added <- 1 / pmax(1, round(1 / added))
    grid <- sort(unique(c(grid, added)))
  }
}

.adaptive_tune <- function(index, ids, candidate, y, d, x, z, response, selection, control, known_response) {
  labels <- if (is.function(candidate$control$penalty_folds)) candidate$control$penalty_folds(ids[index]) else NULL
  plans <- if (selection == "robins") .adaptive_plans(index, ids, labels) else
    list(list(train = index, projection = index, validation = integer(0)))
  designs <- lapply(plans, function(plan) .adaptive_design(z[plan$projection, , drop = FALSE],
    z[if (selection == "robins") plan$validation else plan$train, , drop = FALSE], control))
  cache <- new.env(parent = emptyenv())
  evaluate <- function(scales) vapply(scales, function(scale) {
    key <- format(scale, digits = 17)
    result <- tryCatch({
      pieces <- lapply(seq_along(plans), function(j) {
        plan <- plans[[j]]; scoring_control <- control; scoring_control$seed <- control$seed + j
        yy <- response(plan$train, if (selection == "robins") c(plan$projection, plan$validation) else integer(0))
        model <- .adaptive_fit_at(yy$train, d, x, z, plan$train, candidate, scale, control)
        up <- .adaptive_values(model, x, d, plan$projection)
        # The discrepancy principle is a training-error/noise rule, not CV.
        # It calibrates the actual full-training fit, not a smaller fit whose
        # penalty is later transported to a different sample size. Both its
        # empirical critic and noise envelope use these training rows.
        scored <- if (selection == "robins") plan$validation else plan$train
        uv <- .adaptive_values(model, x, d, scored)
        .adaptive_score(up, uv, if (selection == "robins") yy$validation[seq_along(plan$projection)] else yy$train,
          if (selection == "robins") tail(yy$validation, length(plan$validation)) else yy$train,
          z[plan$projection, , drop = FALSE], z[scored, , drop = FALSE], scoring_control, known_response,
          design = designs[[j]], criterion = selection)
      })
      counts <- vapply(pieces, `[[`, 0L, "n"); factor <- counts / sum(counts)
      score <- sum(factor * vapply(pieces, function(q)
        if (selection == "robins") q$robins[1L, 1L] else q$discrepancy[1L, 1L], 0))
      noise <- Reduce(`+`, Map(function(q, weight) as.vector(q$noise) * weight, pieces, factor))
      list(loss = score, noise = noise, failure = "", projection = lapply(pieces, `[[`, "projection"))
    }, error = function(e) list(loss = Inf, noise = rep(0, control$bootstrap_draws), failure = conditionMessage(e)))
    assign(key, result, cache)
    if (selection == "robins" && !is.finite(result$loss)) stop(result$failure)
    result$loss
  }, 0)
  grid <- .adaptive_grid(candidate, control)
  chosen <- if (selection == "robins") select_penalty_grid(grid, evaluate,
    control$max_extensions, control$extension_factor, control$extension_points) else
      .adaptive_discrepancy_grid(grid, evaluate, cache, control, candidate$method)
  yy <- response(index, integer(0))$train
  fit <- .adaptive_fit_at(yy, d, x, z, index, candidate, chosen$scale, control)
  fit$penalty_cv <- chosen$penalty_cv; fit$penalty_boundary <- chosen$boundary_history
  fit$penalty_status <- chosen$status; fit$selection <- selection
  fit$selection_score_scope <- if (selection == "robins") "held_out_debiased_projection" else "training_discrepancy"
  fit$penalty_parameter <- if (candidate$method == "landweber") "inverse_iterations" else "lambda"
  fit$selection_audit <- plans
  fit$projection_cv <- lapply(ls(cache), function(key) list(scale = as.numeric(key),
    projections = get(key, cache)$projection))
  fit
}

.adaptive_ensemble_weights <- function(G, N, noise, selection, control, solver_control) {
  if (selection == "robins") {
    if (nrow(G) > 10L) stop("The experimental indefinite Robins selector supports at most ten candidates.")
    fit <- .simplex_qp(G, solver_control)
    return(c(fit, list(status = "raw_robins_minimum", threshold = NA_real_, norm_penalty = NA_real_)))
  }
  threshold <- as.numeric(stats::quantile(apply(noise, 1L, max), control$noise_probability, names = FALSE))
  loss <- function(w) as.numeric(crossprod(w, G %*% w))
  minimum <- .simplex_qp(G, solver_control)
  if (loss(minimum$weights) > threshold + 1e-12)
    stop(structure(list(message = "No ensemble mixture meets the calibrated discrepancy threshold.", call = NULL,
      minimum_loss = loss(minimum$weights), threshold = threshold), class = c("cmbridge_discrepancy_error", "error", "condition")))
  limit <- .simplex_qp(N, solver_control)
  if (loss(limit$weights) <= threshold) return(c(limit, list(status = "minimum_norm_feasible",
    threshold = threshold, norm_penalty = Inf)))
  lower <- 0; upper <- max(abs(G)) / max(1e-12, max(abs(N)))
  repeat {
    candidate <- .simplex_qp(G + upper * N, solver_control)
    if (loss(candidate$weights) > threshold) break
    upper <- upper * 10
    if (!is.finite(upper)) stop("Ensemble discrepancy multiplier overflowed.")
  }
  fit <- minimum
  for (iteration in seq_len(60L)) {
    middle <- (lower + upper) / 2
    candidate <- .simplex_qp(G + middle * N, solver_control)
    if (loss(candidate$weights) <= threshold) {lower <- middle; fit <- candidate} else upper <- middle
  }
  c(fit, list(status = "minimum_norm_discrepancy_crossing", threshold = threshold, norm_penalty = lower))
}

.fit_adaptive_ensemble <- function(y, d, x, z, library, fold_id, n_folds, seed,
                                  selection, selection_control, solver_control,
                                  ids = seq_along(d), response = NULL, known_response = FALSE) {
  control <- .adaptive_control(selection_control); control$seed <- seed
  saved_seed <- if (exists(".Random.seed", .GlobalEnv)) get(".Random.seed", .GlobalEnv) else NULL
  on.exit(if (is.null(saved_seed)) {if (exists(".Random.seed", .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)}
    else assign(".Random.seed", saved_seed, .GlobalEnv), add = TRUE)
  set.seed(seed)
  x <- .as_matrix(x); z <- .as_matrix(z); d <- as.numeric(d); y <- as.numeric(y)
  .check_lengths(y, d, x, z)
  if (any(!is.finite(z)) || any(!is.finite(y)) || any(!is.finite(d))) stop("Adaptive input must be finite where observed.")
  if (is.null(response)) response <- function(train, validation) list(train = y[train], validation = y[validation])
  library <- .ensemble_library(library)
  if (any(!vapply(library, function(q) q$method %in% c("sieve_md", "landweber", "pmmr"), TRUE)))
    stop("Adaptive selection supports the three native bridge/adjoint candidates.")
  if (length(ids) != length(y)) stop("selection ids must match input rows.")
  labels <- .ensemble_folds(as.numeric(d != 0), n_folds, fold_id, seed)
  plans <- .adaptive_plans(seq_along(y), ids, labels)
  pieces <- audits <- fold_tuning <- pending <- list(); p <- length(library)
  available <- stats::setNames(rep(TRUE, p), names(library)); failures <- list()
  fit_one <- function(candidate_name, index, scope) {
    if (!available[candidate_name]) return(NULL)
    tryCatch(.adaptive_tune(index, ids, library[[candidate_name]], y, d, x, z,
      response, selection, control, known_response), error = function(e) {
        available[candidate_name] <<- FALSE
        failures[[length(failures) + 1L]] <<- list(candidate = candidate_name,
          scope = scope, training_n = length(index), message = conditionMessage(e), condition = e)
        NULL
      })
  }
  for (j in seq_along(plans)) {
    plan <- plans[[j]]; yy <- response(plan$train, c(plan$projection, plan$validation))
    fits <- stats::setNames(lapply(names(library), fit_one, index = plan$train,
      scope = paste("ensemble fold", j)), names(library))
    up <- matrix(NA_real_, length(plan$projection), p, dimnames = list(NULL, names(library)))
    uv <- matrix(NA_real_, length(plan$validation), p, dimnames = list(NULL, names(library)))
    for (name in names(library)[available]) {
      up[, name] <- .adaptive_values(fits[[name]], x, d, plan$projection)
      uv[, name] <- .adaptive_values(fits[[name]], x, d, plan$validation)
    }
    pending[[j]] <- list(up = up, uv = uv, yy = yy)
    audits[[j]] <- plan; fold_tuning[[j]] <- lapply(fits, `[[`, "penalty_cv")
  }
  final <- stats::setNames(lapply(names(library), fit_one, index = seq_along(y),
    scope = "final training refit"), names(library))
  active_names <- names(library)[available]
  if (!length(active_names)) stop(structure(list(message = "Every adaptive ensemble candidate failed.",
    call = NULL, candidate_failures = failures), class = c("cmbridge_ensemble_candidate_error", "error", "condition")))
  final <- final[active_names]
  for (j in seq_along(plans)) {
    plan <- plans[[j]]; yy <- pending[[j]]$yy
    scoring_control <- control; scoring_control$seed <- control$seed + 1000L + j
    piece <- .adaptive_score(pending[[j]]$up[, active_names, drop = FALSE],
      pending[[j]]$uv[, active_names, drop = FALSE],
      yy$validation[seq_along(plan$projection)], tail(yy$validation, length(plan$validation)),
      z[plan$projection, , drop = FALSE], z[plan$validation, , drop = FALSE], scoring_control, known_response,
      criterion = selection)
    pieces[[j]] <- piece
  }
  counts <- vapply(pieces, `[[`, 0L, "n"); factor <- counts / sum(counts)
  G <- Reduce(`+`, Map(function(q, w) (if (selection == "robins") q$robins else q$discrepancy) * w, pieces, factor))
  N <- Reduce(`+`, Map(function(q, w) q$norm * w, pieces, factor))
  noise <- Reduce(`+`, Map(function(q, w) q$noise * w, pieces, factor))
  solution <- tryCatch(.adaptive_ensemble_weights(G, N, noise, selection, control, solver_control),
    error = function(e) {
      e$candidates <- final; e$candidate_failures <- failures; e$fold_penalty_cv <- fold_tuning
      e$selection_audit <- audits; e$fold_scores <- pieces; e$fold_n_scored <- counts; e$gram <- G
      stop(e)
    })
  active_weights <- stats::setNames(solution$weights, active_names)
  weights <- stats::setNames(rep(0, p), names(library)); weights[active_names] <- active_weights
  used <- active_names[active_weights > 0]
  predict_fun <- .prediction_closure(function(newdata) {
    newdata <- .as_matrix(newdata)
    predictions <- vapply(models, function(model) predict(model, newdata), numeric(nrow(newdata)))
    as.numeric(matrix(predictions, ncol = length(coefficient)) %*% coefficient)
  }, list(models = final[used], coefficient = active_weights[used]))
  fitted <- rep(NA_real_, length(y)); active <- which(d != 0)
  fitted[active] <- predict_fun(x[active, , drop = FALSE])
  scoring_kernel <- .ensemble_kernel(z, "rbf", list(), seed, "u")
  residual <- d * ifelse(is.na(fitted), 0, fitted) - y
  training_loss <- .ensemble_gram(residual, z, scoring_kernel)[1L, 1L]
  dimnames(G) <- list(active_names, active_names)
  out <- list(method = "ensemble", type = "generic", selection = selection, selection_control = control,
    n = length(y), candidates = final, weights = weights, library = library,
    candidate_failures = failures, active_candidates = active_names, fold_id = labels, selection_audit = audits,
    fold_penalty_cv = fold_tuning, fold_scores = pieces, fold_n_scored = counts,
    gram = G, raw_gram = G, norm_gram = N, candidate_cv_loss = diag(G),
    cv_loss = as.numeric(crossprod(active_weights, G %*% active_weights)), solver = solution,
    scoring_kernel = scoring_kernel, moment_loss = training_loss, statistic = "u",
    predict_fun = predict_fun, fitted = fitted, response = y, diagonal = d,
    target_dim = ncol(x), instrument_dim = ncol(z),
    psd_projection = list(applied = FALSE, adjustment_norm = 0),
    source = if (selection == "robins") "https://arxiv.org/abs/2505.20787" else "https://arxiv.org/abs/2603.01337",
    adaptation = "Native objectives retained; estimated conditional projections and bootstrap discrepancy calibration.")
  class(out) <- c("cmbridge_ensemble", "cmbridge_fit")
  out
}
