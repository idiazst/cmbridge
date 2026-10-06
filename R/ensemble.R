# Cross-validated convex stacking for the conditional moments in Section 5.3.2.
# All candidates use the same scoring kernel within each validation fold.

.ensemble_library <- function(library) {
  methods <- c("sieve_md", "landweber", "pmmr", "saturated_l1")
  if (is.null(library)) library <- stats::setNames(as.list(methods[1:3]), methods[1:3])
  if (is.character(library)) library <- as.list(library)
  if (!is.list(library) || !length(library)) {
    stop("library must be a nonempty list of candidate specifications.", call. = FALSE)
  }
  out <- lapply(library, function(x) {
    if (is.character(x) && length(x) == 1L) x <- list(method = x)
    if (!is.list(x) || !is.character(x$method) || length(x$method) != 1L ||
        !x$method %in% methods || length(setdiff(names(x), c("method", "control")))) {
      stop("each candidate must specify method and optional control.", call. = FALSE)
    }
    if (is.null(x$control)) x$control <- list()
    if (!is.list(x$control)) stop("candidate control must be a list.", call. = FALSE)
    x
  })
  nm <- names(out)
  if (is.null(nm)) nm <- rep("", length(out))
  for (j in seq_along(out)) if (is.na(nm[j]) || !nzchar(nm[j])) {
    nm[j] <- paste0(out[[j]]$method, "_", j)
  }
  if (anyDuplicated(nm)) stop("candidate names must be unique.", call. = FALSE)
  names(out) <- nm
  out
}

.ensemble_integer <- function(x, name, minimum = 1L) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) ||
      x != floor(x) || x < minimum || x > .Machine$integer.max) {
    stop(name, " must be an integer >= ", minimum, ".", call. = FALSE)
  }
  as.integer(x)
}

.ensemble_folds <- function(M, n_folds, fold_id, seed) {
  n <- length(M)
  if (!is.null(fold_id)) {
    if (!is.numeric(fold_id) || length(fold_id) != n || anyNA(fold_id) ||
        any(!is.finite(fold_id)) || any(fold_id != floor(fold_id))) {
      stop("fold_id must contain one integer fold label per observation.", call. = FALSE)
    }
    labels <- sort(unique(fold_id))
    if (length(labels) < 2L) stop("at least two folds are required.", call. = FALSE)
    return(match(fold_id, labels))
  }
  k <- .ensemble_integer(n_folds, "n_folds", 2L)
  if (k > n || sum(M == 1) < k) {
    stop("n_folds cannot exceed the number of complete cases.", call. = FALSE)
  }
  set.seed(seed)
  out <- integer(n)
  for (level in c(0, 1)) {
    index <- which(M == level)
    if (length(index)) {
      index <- index[sample.int(length(index))]
      out[index] <- rep(seq_len(k), length.out = length(index))
    }
  }
  out
}

.ensemble_kernel <- function(z, kernel, control, seed) {
  ctrl <- .merge_control(list(
    scale = TRUE, bandwidth = NULL, approximation = "nystrom",
    n_centers = 100L, block_size = 512L, seed = NULL,
    nystrom_ridge = 1e-10, nystrom_tol = 1e-10
  ), control)
  if (!is.logical(ctrl$scale) || length(ctrl$scale) != 1L || is.na(ctrl$scale)) {
    stop("kernel scale must be TRUE or FALSE.", call. = FALSE)
  }
  ctrl$block_size <- .ensemble_integer(ctrl$block_size, "block_size")
  spec <- list(kernel = kernel, control = ctrl)
  spec$scale <- if (ctrl$scale) .scale_fit(z) else
    list(center = rep(0, ncol(z)), scale = rep(1, ncol(z)))
  if (kernel %in% c("linear", "cell")) return(spec)
  zs <- .scale_apply(z, spec$scale)
  kernel_seed <- if (is.null(ctrl$seed)) seed else .ensemble_integer(ctrl$seed, "kernel seed", 0L)
  bw <- ctrl$bandwidth
  if (is.null(bw)) bw <- .median_distance(zs, seed = kernel_seed)
  if (length(bw) != 1L || !is.finite(bw) || bw <= 0) {
    stop("kernel bandwidth must be positive and finite.", call. = FALSE)
  }
  if (length(ctrl$approximation) != 1L || !ctrl$approximation %in% c("nystrom", "exact")) {
    stop("kernel approximation must be 'nystrom' or 'exact'.", call. = FALSE)
  }
  spec$bandwidth <- bw
  if (ctrl$approximation == "nystrom") {
    nc <- .ensemble_integer(ctrl$n_centers, "n_centers")
    centers <- .select_centers(zs, nc, seed = kernel_seed)
    spec$nystrom <- .nystrom_spec(centers, bw, ctrl$nystrom_ridge, ctrl$nystrom_tol)
  }
  spec
}

.ensemble_features <- function(z, spec) {
  zs <- .scale_apply(z, spec$scale)
  if (spec$kernel == "linear") return(cbind(1, zs))
  .nystrom_apply(zs, spec$nystrom)
}

.ensemble_gram <- function(residual, z, spec) {
  residual <- as.matrix(residual)
  n <- nrow(residual)
  if (!n || nrow(z) != n || any(!is.finite(residual))) {
    stop("scoring requires finite residuals on a nonempty validation sample.", call. = FALSE)
  }
  if (spec$kernel == "cell") return(cell_moment_gram(residual, z))
  if (n < 2L) return(matrix(0, ncol(residual), ncol(residual)))
  starts <- seq.int(1L, n, by = spec$control$block_size)
  if (spec$kernel == "linear" || !is.null(spec$nystrom)) {
    moment <- NULL
    diagonal <- matrix(0, ncol(residual), ncol(residual))
    for (start in starts) {
      ii <- seq.int(start, min(n, start + spec$control$block_size - 1L))
      features <- .ensemble_features(z[ii, , drop = FALSE], spec)
      term <- crossprod(features, residual[ii, , drop = FALSE])
      diagonal <- diagonal + crossprod(residual[ii, , drop = FALSE] * sqrt(rowSums(features^2)))
      moment <- if (is.null(moment)) term else moment + term
    }
    G <- (crossprod(moment) - diagonal) / (n * (n - 1))
  } else {
    # Exact equation (22), evaluated in blocks to avoid storing an n-by-n K.
    zs <- .scale_apply(z, spec$scale)
    G <- matrix(0, ncol(residual), ncol(residual))
    for (start in starts) {
      ii <- seq.int(start, min(n, start + spec$control$block_size - 1L))
      K <- rbf_kernel(zs[ii, , drop = FALSE], zs, spec$bandwidth)
      G <- G + crossprod(residual[ii, , drop = FALSE], K %*% residual)
    }
    G <- (G - crossprod(residual)) / (n * (n - 1))
  }
  (G + t(G)) / 2
}

.project_simplex <- function(x) {
  u <- sort(x, decreasing = TRUE)
  j <- seq_along(u)
  rho <- max(which(u - (cumsum(u) - 1) / j > 0))
  out <- pmax(x - (sum(u[seq_len(rho)]) - 1) / rho, 0)
  out / sum(out)
}

.simplex_qp <- function(G, control = list()) {
  ctrl <- .merge_control(list(tol = 1e-8, max_iter = 50000L), control)
  if (length(ctrl$tol) != 1L || !is.finite(ctrl$tol) || ctrl$tol <= 0) {
    stop("solver tol must be positive and finite.", call. = FALSE)
  }
  max_iter <- .ensemble_integer(ctrl$max_iter, "max_iter")
  p <- nrow(G)
  scale <- max(abs(G))
  if (scale == 0) return(list(weights = rep(1 / p, p), kkt_gap = 0, iterations = 0L))
  A <- (G + t(G)) / (2 * scale)
  gap <- function(w) {
    g <- as.numeric(A %*% w)
    max(0, sum(w * g) - min(g))
  }
  w <- rep(0, p)
  w[which.min(diag(A))] <- 1
  if (p <= 10L) {
    # Enumerate faces for small libraries. Singular faces have minimizers on
    # a smaller face, so no ridge or strict positivity constraint is needed.
    best <- as.numeric(crossprod(w, A %*% w))
    for (mask in seq_len(2^p - 1L)) {
      ii <- which(as.logical(intToBits(mask)[seq_len(p)]))
      q <- length(ii)
      lhs <- rbind(cbind(A[ii, ii, drop = FALSE], 1), c(rep(1, q), 0))
      sol <- tryCatch(solve(lhs, c(rep(0, q), 1)), error = function(e) NULL)
      if (is.null(sol) || any(!is.finite(sol)) || any(sol[seq_len(q)] < -ctrl$tol)) next
      candidate <- numeric(p)
      candidate[ii] <- pmax(sol[seq_len(q)], 0)
      candidate <- candidate / sum(candidate)
      value <- as.numeric(crossprod(candidate, A %*% candidate))
      if (value < best) {
        w <- candidate
        best <- value
      }
    }
    if (gap(w) <= ctrl$tol) return(list(weights = w, kkt_gap = gap(w), iterations = 0L))
  }
  L <- max(eigen(A, symmetric = TRUE, only.values = TRUE)$values)
  y <- w
  acceleration <- 1
  for (iter in seq_len(max_iter)) {
    next_w <- .project_simplex(y - as.numeric(A %*% y) / L)
    if (sum(next_w * (A %*% next_w)) > sum(w * (A %*% w)) + 1e-14) {
      y <- w
      acceleration <- 1
      next_w <- .project_simplex(y - as.numeric(A %*% y) / L)
    }
    if (gap(next_w) <= ctrl$tol) {
      return(list(weights = next_w, kkt_gap = gap(next_w), iterations = iter))
    }
    next_acceleration <- (1 + sqrt(1 + 4 * acceleration^2)) / 2
    y <- next_w + (acceleration - 1) / next_acceleration * (next_w - w)
    w <- next_w
    acceleration <- next_acceleration
  }
  stop("simplex optimization did not converge; increase max_iter.", call. = FALSE)
}

.ensemble_loading <- function(phi, train, validation, n) {
  if (is.function(phi)) {
    out <- phi(train, validation)
    if (!is.list(out) || !all(c("train", "validation") %in% names(out)) ||
        length(out$train) != length(train) || length(out$validation) != length(validation)) {
      stop("phi callback must return train and validation vectors matching the supplied indices.", call. = FALSE)
    }
    return(lapply(out[c("train", "validation")], as.numeric))
  }
  if (!is.numeric(phi) || length(phi) != n) {
    stop("phi must be a numeric vector of length n or a loading callback.", call. = FALSE)
  }
  list(train = phi[train], validation = phi[validation])
}

.fit_ensemble <- function(B, V, M, phi, type, library, n_folds, fold_id,
                          seed, kernel, kernel_control, solver_control) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) saved_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) assign(".Random.seed", saved_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  seed <- .ensemble_integer(seed, "seed", 0L)
  B <- .as_matrix(B, "B")
  V <- .as_matrix(V, "V")
  M <- as.numeric(M)
  n <- length(M)
  if (nrow(B) != n || nrow(V) != n || n < 2L || anyNA(M) || any(!M %in% c(0, 1))) {
    stop("B, V, and binary M must have matching lengths, with at least two observations.", call. = FALSE)
  }
  complete <- which(M == 1)
  if (length(complete) < 2L || any(!is.finite(V[complete, , drop = FALSE])) ||
      any(!is.finite(B[if (type == "bridge") seq_len(n) else complete, , drop = FALSE]))) {
    stop("at least two complete cases and finite observed B and V are required.", call. = FALSE)
  }
  library <- .ensemble_library(library)
  p <- length(library)
  folds <- .ensemble_folds(M, n_folds, fold_id, seed)
  k <- max(folds)
  Z <- if (type == "bridge") B else V
  X <- if (type == "bridge") V else B
  cv_predictions <- cv_residuals <- matrix(NA_real_, n, p, dimnames = list(NULL, names(library)))
  grams <- scoring_specs <- scoring_rows <- vector("list", k)
  counts <- integer(k)
  fold_penalty_cv <- fold_penalty_boundary <- vector("list", k)
  available <- stats::setNames(rep(TRUE, p), names(library))
  failures <- list(); fit_scope <- ""
  record_failure <- function(name, message, index, condition = NULL) {
    available[name] <<- FALSE
    failures[[length(failures) + 1L]] <<- list(candidate = name, scope = fit_scope,
      training_n = length(index), message = message, condition = condition)
  }
  fit_candidates <- function(index, loading = NULL) {
    out <- lapply(names(library), function(name) {
      if (!available[name]) return(NULL)
      candidate <- library[[name]]
      tryCatch({
      scales <- candidate$control$penalty_scales
      fit_at <- function(rows, response, control) {
        if (type == "bridge") fit_bridge(B[rows, , drop = FALSE], V[rows, , drop = FALSE],
          M[rows], candidate$method, control)
        else fit_adjoint(B[rows, , drop = FALSE], V[rows, , drop = FALSE],
          M[rows], response, candidate$method, control)
      }
      if (is.null(scales)) return(fit_at(index, loading, candidate$control))
      if (!candidate$method %in% c("sieve_md", "saturated_l1") || !length(scales) ||
          any(!is.finite(scales)) || any(scales <= 0)) {
        stop("Positive penalty_scales are supported for sieve_md and saturated_l1 candidates.")
      }
      ids <- candidate$control$penalty_ids
      if (is.null(ids)) ids <- seq_len(n)
      if (length(ids) != n) stop("penalty_ids must match the ensemble input rows.")
      # Stable person labels agree with the regression/classification learner.
      values <- as.character(ids[index])
      codes <- suppressWarnings(abs(as.numeric(values)))
      nonnumeric <- which(!is.finite(codes))
      codes[nonnumeric] <- vapply(values[nonnumeric], function(id)
        strtoi(substr(digest::digest(id, algo = "xxhash32"), 1, 7), 16L), numeric(1))
      labels <- rep(1L, length(codes))
      # A deeper fit may already contain only one parity group. Use the
      # first binary digit that varies in this training sample, in both
      # packages, without consulting any outcomes or validation data.
      for (bit in 0:52) {
        labels <- as.integer(floor(codes / 2^bit) %% 2) + 1L
        if (length(unique(labels)) == 2L) break
      }
      if (is.function(candidate$control$penalty_folds))
        labels <- candidate$control$penalty_folds(ids[index])
      if (length(labels) != length(index) || anyNA(labels) ||
          any(labels < 1 | labels != as.integer(labels)))
        stop("penalty_folds must return one positive integer label per input row.")
      if (length(unique(labels)) < 2L) stop("Penalty cross-validation needs two person groups.")
      if (any(labels != labels[match(values, values)]))
        stop("A person cannot appear in both penalty training and validation groups.")
      # Construct the loading once per training-only split, then retain it
      # while extending the penalty grid. Conditioning cells always use all Z.
      plans <- lapply(sort(unique(labels)), function(label) {
        train <- index[labels != label]; validation <- index[labels == label]
        observed_validation <- validation[M[validation] == 1]
        if (!any(M[train] == 1) || (type == "adjoint" && !length(observed_validation)))
          stop("Each penalty training and validation group needs measured observations.")
        training_z <- if (type == "bridge") train else train[M[train] == 1]
        list(train = train, validation = validation, observed = observed_validation,
          scoring = .ensemble_kernel(Z[training_z, , drop = FALSE], kernel, kernel_control, seed),
          response = if (type == "adjoint") .ensemble_loading(phi, train, validation, n) else NULL)
      })
      count <- sum(vapply(plans, function(plan)
        length(if (type == "bridge") plan$validation else plan$observed), 0L))
      paired_count <- sum(vapply(plans, function(plan) {
        scored <- if (type == "bridge") plan$validation else plan$observed
        frequencies <- table(.cell_keys(Z[scored, , drop = FALSE]))
        as.integer(sum(frequencies[frequencies > 1L]))
      }, 0L))
      evaluate <- function(requested) vapply(requested, function(scale) {
        loss <- 0
        for (plan in plans) {
          control <- .penalty_control(candidate$control, candidate$method, scale,
                                      length(unique(ids[plan$train])))
          fitted <- fit_at(plan$train, if (type == "adjoint") plan$response$train else NULL, control)
          predictions <- predict(fitted, X[plan$observed, , drop = FALSE])
          if (type == "bridge") {
            residual <- rep(-1, length(plan$validation))
            residual[M[plan$validation] == 1] <- predictions - 1
            scored <- plan$validation
          } else {
            residual <- predictions - plan$response$validation[M[plan$validation] == 1]
            scored <- plan$observed
          }
          loss <- loss + length(scored) * .ensemble_gram(residual, Z[scored, , drop = FALSE], plan$scoring)[1, 1]
        }
        loss / count
      }, 0)
      options <- candidate$control
      selected <- select_penalty_grid(scales, evaluate,
        max_extensions = if (is.null(options$penalty_max_extensions)) 16L else options$penalty_max_extensions,
        extension_factor = if (is.null(options$penalty_extension_factor)) 10 else options$penalty_extension_factor,
        extension_points = if (is.null(options$penalty_extension_points)) 4L else options$penalty_extension_points)
      control <- .penalty_control(candidate$control, candidate$method, selected$scale,
                                  length(unique(ids[index])))
      out <- fit_at(index, loading, control)
      out$penalty_cv <- selected$penalty_cv
      out$penalty_cv$scored_rows <- count
      out$penalty_cv$paired_rows <- paired_count
      out$penalty_cv$informative <- paired_count > 0L
      out$penalty_boundary <- selected$boundary_history
      out$penalty_status <- selected$status
      out$penalty_parameter <- if (candidate$method == "sieve_md") "lambda" else "penalty"
      out$penalty_ids <- ids[index]
      out$penalty_fold_id <- labels
      out
      }, error = function(e) {
        record_failure(name, conditionMessage(e), index, e)
        NULL
      })
    })
    names(out) <- names(library)
    out
  }
  for (fold in seq_len(k)) {
    train <- which(folds != fold)
    validation <- which(folds == fold)
    observed_validation <- validation[M[validation] == 1]
    observed_train <- train[M[train] == 1]
    if (!length(observed_train) || (type == "adjoint" && !length(observed_validation))) {
      stop("each fold must have complete cases in both training and validation samples.", call. = FALSE)
    }
    loading <- if (type == "adjoint") .ensemble_loading(phi, train, validation, n) else NULL
    fit_scope <- paste("validation fold", fold)
    fits <- fit_candidates(train, if (type == "adjoint") loading$train else NULL)
    fold_penalty_cv[[fold]] <- lapply(fits, `[[`, "penalty_cv")
    fold_penalty_boundary[[fold]] <- lapply(fits, `[[`, "penalty_boundary")
    predictions <- matrix(NA_real_, length(observed_validation), p,
      dimnames = list(NULL, names(library)))
    for (name in names(library)[available]) {
      value <- tryCatch(predict(fits[[name]], X[observed_validation, , drop = FALSE]), error = identity)
      if (inherits(value, "error") || any(!is.finite(value))) {
        record_failure(name, if (inherits(value, "error")) conditionMessage(value) else
          "Nonfinite validation predictions.", train, if (inherits(value, "error")) value else NULL)
      } else predictions[, name] <- value
    }
    cv_predictions[observed_validation, ] <- predictions
    if (type == "bridge") {
      residual <- matrix(-1, length(validation), p)
      residual[M[validation] == 1, ] <- predictions - 1
      scored <- validation
    } else {
      phi_validation <- loading$validation[M[validation] == 1]
      if (any(!is.finite(phi_validation))) stop("phi must be finite on validation complete cases.", call. = FALSE)
      residual <- predictions - phi_validation
      scored <- observed_validation
    }
    colnames(residual) <- names(library)
    cv_residuals[scored, ] <- residual
    kernel_train <- if (type == "bridge") train else observed_train
    spec <- .ensemble_kernel(Z[kernel_train, , drop = FALSE], kernel, kernel_control, seed + fold)
    scoring_specs[[fold]] <- spec; scoring_rows[[fold]] <- scored
    # Reject a candidate whose quadratic score overflows, retaining all other
    # candidates on these same validation rows and person assignments.
    for (name in names(library)[available]) {
      score <- tryCatch(.ensemble_gram(residual[, name, drop = FALSE],
        Z[scored, , drop = FALSE], spec), error = identity)
      if (inherits(score, "error") || any(!is.finite(score))) record_failure(name,
        if (inherits(score, "error")) conditionMessage(score) else "Nonfinite validation moment score.",
        train, if (inherits(score, "error")) score else NULL)
    }
    counts[fold] <- length(scored)
  }
  fit_scope <- "final training refit"
  final_loading <- if (type == "adjoint") .ensemble_loading(phi, seq_len(n), integer(0), n)$train else NULL
  final_fits <- fit_candidates(seq_len(n), final_loading)
  active <- names(library)[available]
  if (!length(active)) stop(structure(list(message = "Every conditional-moment ensemble candidate failed.",
    call = NULL, candidate_failures = failures), class = c("cmbridge_ensemble_candidate_error", "error", "condition")))
  final_fits <- final_fits[active]
  for (fold in seq_len(k)) {
    scored <- scoring_rows[[fold]]
    grams[[fold]] <- .ensemble_gram(cv_residuals[scored, active, drop = FALSE],
      Z[scored, , drop = FALSE], scoring_specs[[fold]])
  }
  G <- Reduce(`+`, Map(function(g, count) g * count / sum(counts), grams, counts))
  dimnames(G) <- list(active, active)
  raw_gram <- G
  projection <- .project_gram_psd(raw_gram)
  G <- projection$gram
  solution <- .simplex_qp(G, solver_control)
  active_weights <- stats::setNames(solution$weights, active)
  weights <- stats::setNames(rep(0, p), names(library)); weights[active] <- active_weights
  raw_losses <- projected_losses <- stats::setNames(rep(Inf, p), names(library))
  raw_losses[active] <- diag(raw_gram); projected_losses[active] <- diag(G)
  used <- active[active_weights > 0]
  predict_fun <- .prediction_closure(function(newdata) {
      newdata <- .as_matrix(newdata)
      predictions <- vapply(models, function(fit) predict(fit, newdata), numeric(nrow(newdata)))
      as.numeric(matrix(predictions, nrow = nrow(newdata), ncol = columns) %*% coefficient)
    }, list(models = final_fits[used], coefficient = weights[used], columns = length(used)))
  fitted <- rep(NA_real_, n)
  fitted[complete] <- predict_fun(X[complete, , drop = FALSE])
  training_index <- if (type == "bridge") seq_len(n) else complete
  residual <- if (type == "bridge") M * ifelse(is.na(fitted), 0, fitted) - 1 else
    fitted[complete] - final_loading[complete]
  scoring_kernel <- .ensemble_kernel(Z[training_index, , drop = FALSE], kernel, kernel_control, seed)
  training_loss <- .ensemble_gram(matrix(residual, ncol = 1L), Z[training_index, , drop = FALSE], scoring_kernel)[1, 1]
  out <- list(
    method = "ensemble", type = type, n = n, n_scored = sum(counts),
    library = library, candidates = final_fits, weights = weights,
    active_candidates = active, candidate_failures = failures,
    fold_id = folds, fold_gram = grams, fold_n_scored = counts, gram = G,
    raw_gram = raw_gram, psd_projection = projection[c("eigenvalues",
      "removed_negative_eigenvalues", "adjustment_norm")],
    fold_penalty_cv = fold_penalty_cv, fold_penalty_boundary = fold_penalty_boundary,
    candidate_raw_cv_loss = raw_losses,
    raw_cv_loss = as.numeric(crossprod(active_weights, raw_gram %*% active_weights)),
    cv_predictions = cv_predictions, cv_residuals = cv_residuals,
    candidate_cv_loss = projected_losses,
    cv_loss = as.numeric(crossprod(active_weights, G %*% active_weights)),
    solver = solution[c("kkt_gap", "iterations")], scoring_kernel = scoring_kernel,
    fitted = fitted, moment_loss = training_loss, predict_fun = predict_fun,
    target_dim = ncol(X), instrument_dim = ncol(Z),
    response = if (type == "bridge") rep(1, n) else final_loading[complete],
    diagonal = if (type == "bridge") M else rep(1, length(complete)),
    complete_case_index = complete
  )
  class(out) <- c("cmbridge_ensemble", "cmbridge_fit")
  out
}

#' Cross-validated ensemble for a bridge equation
#' @export
fit_bridge_ensemble <- function(B, V, M, library = NULL, n_folds = 5L,
                                fold_id = NULL, seed = 1L,
                                kernel = c("rbf", "linear", "cell"),
                                kernel_control = list(), solver_control = list()) {
  fit <- .fit_ensemble(B, V, M, NULL, "bridge", library, n_folds, fold_id,
                       seed, match.arg(kernel), kernel_control, solver_control)
  fit$call <- match.call()
  fit
}

#' Cross-validated ensemble for an adjoint equation
#' @export
fit_adjoint_ensemble <- function(B, V, M, phi, library = NULL, n_folds = 5L,
                                 fold_id = NULL, seed = 1L,
                                 kernel = c("rbf", "linear", "cell"),
                                 kernel_control = list(), solver_control = list()) {
  fit <- .fit_ensemble(B, V, M, phi, "adjoint", library, n_folds, fold_id,
                       seed, match.arg(kernel), kernel_control, solver_control)
  fit$call <- match.call()
  fit
}

#' @export
print.cmbridge_ensemble <- function(x, ...) {
  cat("Cross-validated conditional-moment ensemble\n")
  cat("  type: ", x$type, "\n", sep = "")
  cat("  n: ", x$n, "; folds: ", max(x$fold_id), "\n", sep = "")
  cat("  cross-validated moment loss: ", format(x$cv_loss, digits = 5), "\n", sep = "")
  print(x$weights)
  invisible(x)
}
