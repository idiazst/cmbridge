test_that("sieve bridge recovers an in-class polynomial", {
  d <- make_small_bridge()
  f <- fit_bridge(d$B, d$Vobs, d$M, "sieve_md",
                  list(target_basis="poly", instrument_basis="poly",
                       target_degree=3L, instrument_degree=7L, lambda=1e-8))
  g <- seq(-0.9, 0.9, length.out=51)
  expect_lt(sqrt(mean((predict(f,g)-truth_poly(g))^2)), 0.08)
})

test_that("Landweber bridge recovers an in-class polynomial", {
  d <- make_small_bridge()
  f <- fit_bridge(d$B, d$Vobs, d$M, "landweber",
                  list(target_basis="poly", instrument_basis="poly",
                       target_degree=3L, instrument_degree=7L, n_iter=1500L, tol=0))
  g <- seq(-0.9, 0.9, length.out=51)
  expect_lt(sqrt(mean((predict(f,g)-truth_poly(g))^2)), 0.08)
})

test_that("bridge learners accept V missing exactly when M=0", {
  d <- make_small_bridge(n=2000)
  expect_silent(fit_bridge(d$B, d$Vobs, d$M, "sieve_md",
                           list(target_basis="poly", instrument_basis="poly",
                                target_degree=3L, instrument_degree=5L)))
})
