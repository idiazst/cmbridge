test_that("unprojected weights globally minimize an indefinite simplex criterion", {
  G <- matrix(c(-2, .8, 1.4, .8, .2, -.7, 1.4, -.7, -.4), 3L)
  solve_raw <- function(G) cmbridge:::.simplex_qp(G, indefinite = TRUE)$weights
  w <- solve_raw(G)
  grid <- expand.grid(a = seq(0, 1, length.out = 301), b = seq(0, 1, length.out = 301))
  grid <- grid[grid$a + grid$b <= 1, ]
  W <- cbind(grid$a, grid$b, 1 - grid$a - grid$b)
  expect_lte(drop(crossprod(w, G %*% w)), min(rowSums((W %*% G) * W)) + 1e-12)
  expect_equal(sum(w), 1)
  expect_true(all(w >= 0))
  # A constant offset to all entries leaves the feasible objective ordering unchanged.
  expect_equal(solve_raw(G + 20), w, tolerance = 1e-10)
  expect_error(cmbridge:::.simplex_qp(diag(11), indefinite = TRUE), "at most 10")
})

test_that("PSD switch preserves fitted candidates, penalty CV and raw scores", {
  d <- make_small_bridge(n = 300, seed = 81)
  folds <- rep(1:3, length.out = length(d$M))
  lib <- list(linear = list(method = "sieve_md", control = list(
    target_basis = "poly", instrument_basis = "poly", target_degree = 1L,
    instrument_degree = 2L, penalty_scales = 10^seq(-8, 0, by = 2))),
    quadratic = list(method = "sieve_md", control = list(
      target_basis = "poly", instrument_basis = "poly", target_degree = 2L,
      instrument_degree = 3L, lambda = .01)))
  for (type in c("bridge", "adjoint")) {
    fit <- function(project) if (type == "bridge") fit_bridge_ensemble(
      d$B, d$Vobs, d$M, lib, fold_id = folds, project_psd = project) else
        fit_adjoint_ensemble(d$B, d$Vobs, d$M, 1 + d$Vobs, lib,
          fold_id = folds, project_psd = project)
    yes <- fit(TRUE); no <- fit(FALSE)
    expect_equal(yes$raw_gram, no$raw_gram, tolerance = 0)
    expect_equal(yes$cv_predictions, no$cv_predictions, tolerance = 0)
    expect_equal(yes$fold_penalty_cv, no$fold_penalty_cv, tolerance = 0)
    expect_equal(lapply(yes$candidates, `[[`, "penalty_cv"),
      lapply(no$candidates, `[[`, "penalty_cv"), tolerance = 0)
    expect_true(yes$psd_projection$applied)
    expect_false(no$psd_projection$applied)
    expect_equal(no$gram, no$raw_gram, tolerance = 0)
    expect_lte(no$raw_cv_loss, yes$raw_cv_loss + 1e-12)
    z <- matrix(seq(-1, 1, length.out = 13))
    for (candidate in names(no$candidates))
      expect_equal(predict(yes$candidates[[candidate]], z),
        predict(no$candidates[[candidate]], z), tolerance = 0)
    expect_equal(predict(no, z), as.numeric(vapply(no$candidates, predict,
      numeric(nrow(z)), newdata = z) %*% no$weights[names(no$candidates)]))
  }
  expect_error(fit_bridge_ensemble(d$B, d$Vobs, d$M, lib, project_psd = NA), "TRUE or FALSE")
})
