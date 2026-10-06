test_that("sieve ridge CV uses raw held-out U-statistics and the supplied person splits", {
  set.seed(947)
  n <- 900L
  H <- rbinom(n, 1, .5); A <- rbinom(n, 1, .5)
  latent <- rbinom(n, 1, .15 + .7*A)
  M <- rbinom(n, 1, .3 + .3*latent)
  B <- cbind(H, A); V <- cbind(H, latent); V[M == 0, ] <- NA_real_
  outer <- rep(1:3, length.out = n)
  calls <- list()
  shared <- function(ids) {
    calls[[length(calls) + 1L]] <<- ids
    (floor((ids - 1)/3) %% 2) + 1L
  }
  library <- list(sieve = list(method = "sieve_md", control = list(
    target_basis = "cell", instrument_basis = "cell", link = "inverse_logit",
    penalty_scales = 10^seq(-10, 0, by = .5), penalty_ids = seq_len(n),
    penalty_folds = shared, penalty_max_extensions = 8L,
    penalty_extension_factor = 10, penalty_extension_points = 4L)))
  f <- fit_bridge_ensemble(B, V, M, library, fold_id = outer, kernel = "cell")
  fit <- f$candidates$sieve; cv <- fit$penalty_cv
  expect_identical(fit$penalty_parameter, "lambda")
  expect_equal(fit$tuning$lambda, cv$scale[cv$selected])
  expect_gt(fit$tuning$lambda, 1e-8)
  expect_true(which(cv$selected) > 1 && which(cv$selected) < nrow(cv))
  expect_true(all(cv$scale > 0))
  expect_length(calls, 4L)
  for (fold in 1:3) expect_identical(calls[[fold]], which(outer != fold))
  expect_identical(calls[[4]], seq_len(n))
  labels <- (floor((seq_len(n) - 1)/3) %% 2) + 1L
  expect_equal(fit$penalty_fold_id, labels)
  expect_identical(fit$penalty_ids, seq_len(n))
  independent <- vapply(cv$scale, function(lambda) {
    loss <- 0
    for (label in 1:2) {
      train <- which(labels != label); valid <- which(labels == label)
      observed <- valid[M[valid] == 1]
      model <- fit_bridge(B[train, ], V[train, ], M[train], "sieve_md",
        list(target_basis = "cell", instrument_basis = "cell", link = "inverse_logit", lambda = lambda))
      residual <- rep(-1, length(valid))
      residual[M[valid] == 1] <- predict(model, V[observed, ]) - 1
      loss <- loss + length(valid) * as.numeric(cell_moment_gram(residual, B[valid, ]))
    }
    loss/n
  }, 0)
  expect_equal(cv$loss, independent, tolerance = 1e-12)
  expect_true(all(is.finite(fitted(f)[M == 1])))

  library$sieve$control$link <- "identity"
  loading_calls <- list()
  phi <- function(train, validation) {
    expect_length(intersect(train, validation), 0L)
    loading_calls[[length(loading_calls) + 1L]] <<- list(train = train, validation = validation)
    list(train = rep(2, length(train)), validation = rep(2, length(validation)))
  }
  adjoint <- fit_adjoint_ensemble(B, V, M, phi, library, fold_id = outer, kernel = "cell")
  a <- adjoint$candidates$sieve
  expect_equal(a$tuning$lambda, a$penalty_cv$scale[a$penalty_cv$selected])
  expect_identical(a$tuning$link, "identity")
  expect_length(loading_calls, 12L)
  expect_true(any(vapply(loading_calls, function(x) length(x$train) < n/2, FALSE)))
  expect_true(all(a$penalty_boundary$status == "boundary_failure"))
  expect_true(which(a$penalty_cv$selected) > 1 && which(a$penalty_cv$selected) < nrow(a$penalty_cv))
})

test_that("ridge and L1 CV controls preserve their different penalty conventions", {
  control <- list(penalty_scales = .01, penalty_ids = 1:10, penalty_folds = identity,
    penalty_max_extensions = 3, penalty_extension_factor = 10, penalty_extension_points = 4,
    target_basis = "cell")
  ridge <- cmbridge:::.penalty_control(control, "sieve_md", .01, 100)
  l1 <- cmbridge:::.penalty_control(control, "saturated_l1", .01, 100)
  expect_equal(ridge, list(target_basis = "cell", lambda = .01))
  expect_equal(l1, list(target_basis = "cell", penalty = .001))
})
