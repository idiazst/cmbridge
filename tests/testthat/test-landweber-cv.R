test_that("Landweber weighting CV reproduces held-out losses on shared training splits", {
  set.seed(219)
  n <- 720L
  H <- rbinom(n, 1, .5); A <- rbinom(n, 1, .5)
  latent <- rbinom(n, 1, .2 + .6 * A)
  M <- rbinom(n, 1, .35 + .3 * latent)
  B <- cbind(H, A); V <- cbind(H, latent); V[M == 0, ] <- NA_real_
  folds <- rep(1:3, length.out = n)
  calls <- list()
  shared <- function(ids) {
    calls[[length(calls) + 1L]] <<- ids
    (floor((ids - 1) / 3) %% 2) + 1L
  }
  base <- list(target_basis = "poly", target_degree = 1L,
    instrument_basis = "cell", link = "inverse_logit", precondition = TRUE)
  control <- c(base, list(penalty_scales = 10^seq(-4, 0, by = .5),
    penalty_ids = seq_len(n), penalty_folds = shared))
  fit <- fit_bridge_ensemble(B, V, M,
    list(iterative = list(method = "landweber", control = control)),
    fold_id = folds, kernel = "cell")
  candidate <- fit$candidates$iterative; cv <- candidate$penalty_cv
  expect_length(fit$candidate_failures, 0L)
  expect_identical(candidate$penalty_parameter, "weight_ridge")
  expect_equal(candidate$tuning$weight_ridge, cv$scale[cv$selected])
  expect_true(all(cv$scale > 0))
  chosen <- which(cv$selected)
  expect_true(length(chosen) == 1L && chosen > 1L && chosen < nrow(cv))
  expect_equal(cv$loss[chosen], min(cv$loss))
  expect_true(candidate$solver$stopped_by_tolerance)
  expect_length(calls, 4L)
  for (fold in 1:3) expect_identical(calls[[fold]], which(folds != fold))
  expect_identical(calls[[4L]], seq_len(n))
  expect_identical(candidate$penalty_ids, seq_len(n))
  labels <- (floor((seq_len(n) - 1) / 3) %% 2) + 1L
  expect_equal(candidate$penalty_fold_id, labels)
  for (row in which(is.finite(cv$loss))) {
    loss <- 0
    for (label in 1:2) {
      train <- which(labels != label); valid <- which(labels == label)
      settings <- c(base, list(weight_ridge = cv$scale[row]))
      model <- fit_bridge(B[train, ], V[train, ], M[train], "landweber", settings)
      residual <- rep(-1, length(valid)); measured <- M[valid] == 1
      residual[measured] <- predict(model, V[valid[measured], ]) - 1
      # Independent ordered-pair formula, not the scorer used in fitting.
      score <- 0
      for (group in split(seq_along(valid), cmbridge:::.cell_keys(B[valid, ]))) {
        k <- length(group)
        if (k > 1L) score <- score +
          (sum(residual[group])^2 - sum(residual[group]^2)) / (length(valid) * (k - 1L))
      }
      loss <- loss + length(valid) * score
    }
    expect_equal(cv$loss[row], loss / n, tolerance = 1e-12)
  }
  expect_true(all(fitted(candidate)[M == 1] >= 1))
})

test_that("conditioning-weight tuning preserves other control values", {
  control <- list(target_basis = "poly", instrument_basis = "cell", n_iter = 17L,
    penalty_scales = c(.01, .1), penalty_ids = 1:10, penalty_folds = identity)
  tuned <- cmbridge:::.penalty_control(control, "landweber", .1, 100)
  expect_equal(tuned, list(target_basis = "poly", instrument_basis = "cell",
    n_iter = 17L, weight_ridge = .1))
})
