test_that("the common kernel Gram implements equation (22)", {
  set.seed(31)
  z <- matrix(rnorm(60), ncol = 2)
  residual <- matrix(rnorm(90), ncol = 3)
  spec <- cmbridge:::.ensemble_kernel(z, "rbf",
    list(scale = FALSE, bandwidth = 0.7, approximation = "exact", block_size = 7L), 1L)
  direct <- crossprod(residual, rbf_kernel(z, z, 0.7) %*% residual) / nrow(z)^2
  expect_equal(cmbridge:::.ensemble_gram(residual, z, spec), direct, tolerance = 1e-12)
  spec2 <- cmbridge:::.ensemble_kernel(z, "rbf",
    list(scale = FALSE, bandwidth = 0.7, n_centers = nrow(z)), 1L)
  expect_equal(cmbridge:::.ensemble_gram(residual, z, spec2), direct, tolerance = 1e-7)
})

test_that("simplex optimization handles interior, boundary, and singular optima", {
  sol <- cmbridge:::.simplex_qp(diag(c(1, 2, 4)))
  expect_equal(sol$weights, c(4, 2, 1) / 7, tolerance = 1e-8)
  G <- tcrossprod(c(0, 1, 2))
  sol <- cmbridge:::.simplex_qp(G)
  expect_equal(sol$weights, c(1, 0, 0))
  sol <- cmbridge:::.simplex_qp(matrix(1, 3, 3))
  expect_equal(sum(sol$weights), 1)
  expect_equal(sol$kkt_gap, 0)
  expect_equal(cmbridge:::.simplex_qp(matrix(0, 2, 2))$weights, c(0.5, 0.5))
  # More than ten candidates exercises the projected-gradient solver.
  sol <- cmbridge:::.simplex_qp(diag(seq_len(12)))
  expect_equal(sol$weights, (1 / seq_len(12)) / sum(1 / seq_len(12)), tolerance = 1e-6)
})

ensemble_test_library <- function() list(
  linear = list(method = "sieve_md", control = list(target_basis = "poly",
    instrument_basis = "poly", target_degree = 1L, instrument_degree = 5L)),
  cubic = list(method = "sieve_md", control = list(target_basis = "poly",
    instrument_basis = "poly", target_degree = 3L, instrument_degree = 5L))
)

test_that("inner predictions exclude their validation fold and ignore missing V", {
  d <- make_small_bridge(n = 1200)
  folds <- rep(1:3, length.out = length(d$M))
  lib <- ensemble_test_library()
  f <- fit_bridge_ensemble(d$B, d$Vobs, d$M, lib, fold_id = folds,
    kernel_control = list(approximation = "exact", block_size = 64L))
  train <- which(folds != 1)
  validation <- which(folds == 1 & d$M == 1)
  manual <- fit_bridge(d$B[train], d$Vobs[train], d$M[train],
    lib$cubic$method, lib$cubic$control)
  expect_equal(f$cv_predictions[validation, "cubic"], predict(manual, d$Vobs[validation]))
  expect_true(all(is.na(f$cv_predictions[d$M == 0, ])))
  expect_true(all(f$cv_residuals[d$M == 0, ] == -1))
  expect_true(all(f$weights >= 0))
  expect_equal(sum(f$weights), 1)
  expect_lte(f$cv_loss, min(f$candidate_cv_loss) + 1e-12)
  g <- seq(-0.8, 0.8, length.out = 41)
  combined <- Reduce(`+`, Map(function(x, w) predict(x, g) * w, f$candidates, f$weights))
  expect_equal(predict(f, g), combined)
  kernel <- f$scoring_kernel
  z <- cmbridge:::.scale_apply(matrix(d$B, ncol = 1), kernel$scale)
  pred <- rep(0, length(d$M))
  pred[d$M == 1] <- predict(f, d$Vobs[d$M == 1])
  r <- 1 - d$M * pred
  direct_loss <- as.numeric(crossprod(r, rbf_kernel(z, z, kernel$bandwidth) %*% r)) / length(r)^2
  expect_equal(moment_loss(f, rep(1, length(r)), d$M, d$Vobs, d$B), direct_loss)
})

test_that("adjoint loading is common, trained within folds, and scored on complete cases", {
  d <- make_small_bridge(n = 1200, seed = 92)
  B <- d$B
  B[d$M == 0] <- NA_real_
  calls <- list()
  phi <- function(train, validation) {
    calls[[length(calls) + 1L]] <<- list(train = train, validation = validation)
    list(train = 1 + d$Vobs[train], validation = 1 + d$Vobs[validation])
  }
  f <- fit_adjoint_ensemble(B, d$Vobs, d$M, phi, ensemble_test_library(), n_folds = 3L)
  expect_length(calls, 4L)
  for (i in 1:3) expect_length(intersect(calls[[i]]$train, calls[[i]]$validation), 0L)
  expect_equal(calls[[4]]$train, seq_along(d$M))
  expect_length(calls[[4]]$validation, 0L)
  expect_equal(f$n_scored, sum(d$M))
  expect_equal(f$fold_n_scored, as.integer(table(f$fold_id[d$M == 1])))
  expect_true(all(is.na(f$cv_residuals[d$M == 0, ])))
  expect_equal(f$cv_residuals[d$M == 1, ],
    f$cv_predictions[d$M == 1, ] - (1 + d$Vobs[d$M == 1]))
  weights <- f$fold_n_scored / sum(f$fold_n_scored)
  expected <- Reduce(`+`, Map(function(G, w) G * w, f$fold_gram, weights))
  dimnames(expected) <- dimnames(f$gram)
  expect_equal(f$gram, expected)
})

test_that("ensembles validate inputs and preserve the calling random state", {
  d <- make_small_bridge(n = 150)
  set.seed(325)
  before <- .Random.seed
  lib <- ensemble_test_library()
  lib$pmmr <- list(method = "pmmr")
  expect_silent(fit_bridge_ensemble(d$B, d$Vobs, d$M, library = lib, n_folds = 3L))
  expect_identical(.Random.seed, before)
  expect_error(fit_bridge_ensemble(d$B, d$Vobs, d$M, library = list()), "nonempty")
  expect_error(fit_bridge_ensemble(d$B, d$Vobs, rep(0.5, 150)), "binary")
  expect_error(fit_bridge_ensemble(d$B, d$Vobs, d$M, fold_id = rep(1, 150)), "two folds")
  bad_folds <- ifelse(d$M == 1, 1, 2)
  expect_error(fit_bridge_ensemble(d$B, d$Vobs, d$M, fold_id = bad_folds), "complete cases")
  expect_error(fit_adjoint_ensemble(d$B, d$Vobs, d$M,
    function(train, validation) 1, n_folds = 3L), "callback")
})
