test_that("adaptive splits keep fitting, projection and scoring people separate", {
  ids <- rep(seq_len(60), each = 2)
  for (k in 2:3) {
    plans <- cmbridge:::.adaptive_plans(seq_along(ids), ids, rep(rep(1:k, length.out = 60), each = 2))
    for (plan in plans) {
      expect_length(intersect(ids[plan$train], ids[plan$projection]), 0)
      expect_length(intersect(ids[plan$train], ids[plan$validation]), 0)
      expect_length(intersect(ids[plan$projection], ids[plan$validation]), 0)
    }
  }
  expect_equal(cmbridge:::.adaptive_labels(1:120)[c(2, 21, 40)],
    cmbridge:::.adaptive_labels(c(2, 21, 40)))
})

test_that("Robins score is the paper score up to the common constant", {
  set.seed(1)
  zp <- matrix(rep(0:1, 100), ncol = 1)
  zv <- matrix(rep(0:1, 70), ncol = 1)
  up <- cbind(1 + .4 * zp[, 1], 1 - .1 * zp[, 1])
  uv <- cbind(1 + .4 * zv[, 1], 1 - .1 * zv[, 1])
  ctrl <- cmbridge:::.adaptive_control(list(bootstrap_draws = 19L))
  got <- cmbridge:::.adaptive_score(up, uv, rep(1, nrow(up)), rep(1, nrow(uv)), zp, zv, ctrl, TRUE)
  expect_equal(diag(got$robins), colMeans((uv - 1)^2), tolerance = 1e-7)
  expect_true(all(got$projection$ridge > 0))
  expect_false(got$projection$ridge_boundary)
  w <- c(.3, .7)
  combined <- cmbridge:::.adaptive_score(up %*% w, uv %*% w, rep(1, nrow(up)),
    rep(1, nrow(uv)), zp, zv, ctrl, TRUE)
  expect_equal(as.numeric(crossprod(w, got$robins %*% w)), combined$robins[1, 1], tolerance = 1e-7)
  expect_equal(got$noise, cmbridge:::.adaptive_score(up, uv, rep(1, nrow(up)),
    rep(1, nrow(uv)), zp, zv, ctrl, TRUE)$noise)
})

test_that("discrepancy ensemble finds a feasible minimum-norm mixture", {
  G <- diag(c(.01, 1)); N <- diag(c(1, .01)); noise <- matrix(.1, 19, 2)
  got <- cmbridge:::.adaptive_ensemble_weights(G, N, noise, "syrgkanis",
    cmbridge:::.adaptive_control(list(bootstrap_draws = 19L)), list())
  expect_equal(sum(got$weights), 1)
  expect_true(all(got$weights >= 0))
  expect_lte(as.numeric(crossprod(got$weights, G %*% got$weights)), got$threshold + 1e-9)
  expect_error(cmbridge:::.adaptive_ensemble_weights(diag(2), N, noise, "syrgkanis",
    cmbridge:::.adaptive_control(list(bootstrap_draws = 19L)), list()), "No ensemble mixture")
})

test_that("both selectors tune all three bridge and adjoint methods", {
  set.seed(91)
  n <- 600
  B <- matrix(rbinom(n, 1, .5), ncol = 1)
  V <- matrix(rbinom(n, 1, .1 + .8 * B[, 1]), ncol = 1)
  M <- rbinom(n, 1, .5 + .2 * V[, 1])
  missing <- V; missing[M == 0] <- NA_real_
  for (selection in c("robins", "syrgkanis")) for (method in c("sieve_md", "landweber", "pmmr")) {
    native <- list(link = "inverse_logit", penalty_scales = c(1e-5, 1e-3, .1, 1),
      iteration_grid = c(2L, 10L, 100L, 1000L))
    if (method == "pmmr") {
      native$target_centers <- 4L; native$instrument_centers <- 4L
      native$target_centers_matrix <- native$instrument_centers_matrix <- matrix(0:1, ncol = 1)
    } else {
      native$target_basis <- native$instrument_basis <- "poly"
      native$target_degree <- native$instrument_degree <- 1L
    }
    bridge <- fit_bridge(B, missing, M, method, native, selection,
      list(bootstrap_draws = 19L))
    expect_true(all(predict(bridge, V) >= 1))
    expect_equal(bridge$selection, selection)
    expect_true(any(bridge$penalty_cv$selected))
    native$link <- "identity"
    adjoint <- fit_adjoint(B, V, M, .5 + V[, 1], method, native, selection,
      list(bootstrap_draws = 19L))
    expect_true(all(is.finite(predict(adjoint, B))))
    expect_equal(adjoint$selection, selection)
  }
})

test_that("Robins ensembles preserve callbacks, raw scores and vector prediction", {
  set.seed(17); n <- 360
  B <- matrix(rbinom(n,1,.5),ncol=1)
  V <- matrix(rbinom(n,1,.1+.8*B[,1]),ncol=1)
  M <- rbinom(n,1,.6+.1*V[,1]); ids <- seq_len(n)
  library <- list(sieve=list(method="sieve_md",control=list(link="inverse_logit",
    target_basis="poly",instrument_basis="poly",target_degree=1L,instrument_degree=1L,
    penalty_scales=c(1e-5,1e-3,.1,1))))
  bridge <- fit_bridge_ensemble(B,V,M,library,fold_id=ids%%3+1,
    selection="robins",selection_control=list(bootstrap_draws=19L))
  expect_true(all(predict(bridge,0:1)>=1))
  expect_false(bridge$psd_projection$applied)
  expect_equal(bridge$gram,bridge$raw_gram)
  expect_true(is.finite(moment_loss(bridge)))
  expect_true(is.finite(moment_loss(bridge,rep(1,n),M,V,B)))
  callback <- function(train,validation) {
    expect_length(intersect(train,validation),0)
    list(train=(.5+V[,1])[train],validation=(.5+V[,1])[validation])
  }
  library$sieve$control$link <- "identity"
  adjoint <- fit_adjoint_ensemble(B,V,M,callback,library,fold_id=ids%%3+1,
    selection="robins",selection_control=list(bootstrap_draws=19L))
  expect_true(all(is.finite(predict(adjoint,0:1))))
})
