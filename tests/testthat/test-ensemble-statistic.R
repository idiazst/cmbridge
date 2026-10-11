test_that("V kernel scores retain self-products and use the squared sample size", {
  set.seed(913)
  z <- matrix(rnorm(24), ncol = 2L)
  residual <- matrix(rnorm(36), ncol = 3L)
  for (kernel in c("linear", "rbf")) for (approximation in c("exact", "nystrom")) {
    spec <- cmbridge:::.ensemble_kernel(z, kernel, list(scale = FALSE,
      bandwidth = .8, approximation = approximation, n_centers = nrow(z), block_size = 5L),
      1L, statistic = "v")
    K <- if (kernel == "linear") tcrossprod(cbind(1, z)) else rbf_kernel(z, z, .8)
    expected <- crossprod(residual, K %*% residual) / nrow(z)^2
    G <- cmbridge:::.ensemble_gram(residual, z, spec)
    expect_equal(G, expected, tolerance = 1e-7)
    expect_gte(min(eigen(G, symmetric = TRUE, only.values = TRUE)$values), -1e-12)
    expect_equal(cmbridge:::.ensemble_gram(residual[1, , drop = FALSE],
      z[1, , drop = FALSE], spec),
      tcrossprod(residual[1, ]) * K[1, 1], tolerance = 1e-7)
  }
  # With a low-rank approximation, remove its own feature diagonal, not 1.
  spec <- cmbridge:::.ensemble_kernel(z, "rbf", list(bandwidth = .8, n_centers = 3L), 1L)
  F <- cmbridge:::.ensemble_features(z, spec)
  GU <- cmbridge:::.ensemble_gram(residual, z, spec)
  GV <- cmbridge:::.ensemble_gram(residual, z, spec, statistic = "v")
  diagonal <- crossprod(residual * sqrt(rowSums(F^2)))
  expect_equal(GV, (nrow(z) - 1) / nrow(z) * GU + diagonal / nrow(z)^2)
})

test_that("cell V scores are the frequency-weighted squared cell means", {
  z <- matrix(c(0, 0, 1, 1, 1, 2))
  residual <- cbind(c(-1, 3, 2, -2, 1, 4), c(2, 1, 0, 3, -1, 2))
  expected <- Reduce(`+`, lapply(split(seq_len(nrow(z)), z[, 1]), function(ii)
    tcrossprod(colMeans(residual[ii, , drop = FALSE])) * length(ii) / nrow(z)))
  expect_equal(as.numeric(cell_moment_gram(residual, z, "v")), as.numeric(expected))
  singleton <- cell_moment_gram(residual, matrix(seq_len(nrow(z))), "v")
  expect_equal(as.numeric(singleton), as.numeric(crossprod(residual) / nrow(z)))
  expect_equal(as.numeric(cell_moment_gram(residual[1, , drop = FALSE], z[1, , drop = FALSE], "v")),
    as.numeric(tcrossprod(residual[1, ])))
  expect_error(cell_moment_gram(residual, z, "invalid"), "arg")
})

test_that("V bridge and adjoint ensembles keep their PSD scores without projection", {
  d <- make_small_bridge(n = 300, seed = 76)
  folds <- rep(1:3, length.out = length(d$M))
  lib <- list(sieve = list(method = "sieve_md", control = list(
    target_basis = "poly", instrument_basis = "poly", target_degree = 2L, instrument_degree = 3L)))
  for (type in c("bridge", "adjoint")) {
    f <- if (type == "bridge") fit_bridge_ensemble(d$B, d$Vobs, d$M, lib,
      fold_id = folds, statistic = "v") else
      fit_adjoint_ensemble(d$B, d$Vobs, d$M, 1 + d$Vobs, lib, fold_id = folds, statistic = "v")
    expect_identical(f$statistic, "v")
    expect_identical(f$scoring_kernel$statistic, "v")
    expect_false(f$psd_projection$applied)
    expect_equal(f$gram, f$raw_gram)
    expect_equal(f$psd_projection$adjustment_norm, 0)
    expect_gte(f$cv_loss, -1e-12)
    expect_equal(sum(f$weights), 1)
    Z <- matrix(if (type == "bridge") d$B else d$Vobs, ncol = 1L)
    for (fold in 1:3) {
      train <- which(folds != fold)
      vi <- which(folds == fold)
      if (type == "adjoint") {train <- train[d$M[train] == 1]; vi <- vi[d$M[vi] == 1]}
      spec <- cmbridge:::.ensemble_kernel(Z[train, , drop = FALSE], "rbf", list(), 1L + fold)
      F <- cmbridge:::.ensemble_features(Z[vi, , drop = FALSE], spec)
      expected <- crossprod(crossprod(F, f$cv_residuals[vi, , drop = FALSE])) / length(vi)^2
      expect_equal(f$fold_gram[[fold]], expected, ignore_attr = TRUE)
    }
    rows <- if (type == "bridge") seq_along(d$M) else which(d$M == 1)
    y <- if (type == "bridge") rep(1, length(rows)) else 1 + d$Vobs[rows]
    diagonal <- if (type == "bridge") d$M else rep(1, length(rows))
    target <- matrix(if (type == "bridge") d$Vobs[rows] else d$B[rows], ncol = 1L)
    residual <- y
    observed <- diagonal != 0
    residual[observed] <- residual[observed] - diagonal[observed] * predict(f, target[observed, , drop = FALSE])
    F <- cmbridge:::.ensemble_features(Z[rows, , drop = FALSE], f$scoring_kernel)
    expected <- sum(crossprod(F, residual)^2) / length(rows)^2
    expect_equal(moment_loss(f), expected, tolerance = 1e-12)
    expect_equal(moment_loss(f, y, diagonal, target, Z[rows, , drop = FALSE]), expected, tolerance = 1e-12)
  }
  default <- fit_bridge_ensemble(d$B, d$Vobs, d$M, lib, fold_id = folds)
  explicit <- fit_bridge_ensemble(d$B, d$Vobs, d$M, lib, fold_id = folds, statistic = "u")
  expect_equal(default$weights, explicit$weights)
  expect_equal(default$raw_gram, explicit$raw_gram)
  expect_true(default$psd_projection$applied)
  expect_error(fit_bridge_ensemble(d$B, d$Vobs, d$M, lib, statistic = "invalid"), "arg")
})

test_that("V penalty CV uses training-only V scores on the supplied person splits", {
  set.seed(947)
  n <- 600L
  H <- rbinom(n, 1, .5); A <- rbinom(n, 1, .5)
  latent <- rbinom(n, 1, .15 + .7 * A)
  M <- rbinom(n, 1, .3 + .3 * latent)
  B <- cbind(H, A); V <- cbind(H, latent); V[M == 0, ] <- NA_real_
  labels <- rep(1:3, length.out = n)
  penalty_folds <- function(ids) (floor((ids - 1) / 3) %% 2) + 1L
  controls <- list(target_basis = "cell", instrument_basis = "cell", link = "identity",
    penalty_scales = 10^seq(-8, 0, by = 1), penalty_ids = seq_len(n),
    penalty_folds = penalty_folds)
  for (type in c("bridge", "adjoint")) {
    lib <- list(sieve = list(method = "sieve_md", control = controls))
    phi <- 1 + latent
    f <- if (type == "bridge") fit_bridge_ensemble(B, V, M, lib,
      fold_id = labels, kernel = "cell", statistic = "v") else
      fit_adjoint_ensemble(B, V, M, phi, lib, fold_id = labels, kernel = "cell", statistic = "v")
    cv <- f$candidates$sieve$penalty_cv
    expect_true(all(cv$statistic == "v"))
    expect_true(grepl("^interior", f$candidates$sieve$penalty_status))
    folds <- penalty_folds(seq_len(n))
    independent <- vapply(cv$scale, function(lambda) {
      total <- count <- 0
      for (fold in 1:2) {
        tr <- which(folds != fold); vi <- which(folds == fold)
        measured <- vi[M[vi] == 1]
        ctrl <- list(target_basis = "cell", instrument_basis = "cell", link = "identity", lambda = lambda)
        fit <- if (type == "bridge") fit_bridge(B[tr, ], V[tr, ], M[tr], "sieve_md", ctrl) else
          fit_adjoint(B[tr, ], V[tr, ], M[tr], phi[tr], "sieve_md", ctrl)
        if (type == "bridge") {
          r <- rep(-1, length(vi)); r[M[vi] == 1] <- predict(fit, V[measured, ]) - 1
          z <- B[vi, , drop = FALSE]
        } else {r <- predict(fit, B[measured, ]) - phi[measured]; z <- V[measured, , drop = FALSE]}
        total <- total + length(r) * as.numeric(cell_moment_gram(r, z, "v"))
        count <- count + length(r)
      }
      total / count
    }, 0)
    expect_equal(cv$loss, independent, tolerance = 1e-12)
  }
})
