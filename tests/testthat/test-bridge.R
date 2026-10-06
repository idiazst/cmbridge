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

test_that("raw conditional-moment predictions are unbounded and finite bounds fail visibly", {
  x <- matrix(seq(-1, 1, length.out = 100), ncol = 1)
  y <- 4 * x[, 1]
  for (method in c("sieve_md", "landweber", "pmmr", "saturated_l1")) {
    f <- fit_cm(y, target = x, instrument = x, method = method)
    expect_lt(min(predict(f, x)), 0)
    expect_equal(f$residual, y - fitted(f))
    expect_error(fit_cm(y, target = x, instrument = x, method = method,
      control = list(lower = 1, upper = 2)), "Post-fit prediction bounds")
    expect_error(fit_cm(y, target = x, instrument = x, method = method,
      control = list(lower = 2, upper = 1)), "ordered scalar")
  }
})

test_that("joint-category sieve represents interactions and has an unseen-cell fallback", {
  grid <- expand.grid(x1 = 0:1, x2 = 0:1)
  x <- as.matrix(grid[rep(seq_len(nrow(grid)), each = 20), ])
  y <- 1 + x[, 1] * x[, 2] + 2 * x[, 2]
  f <- fit_cm(y, target = x, instrument = x, method = "sieve_md",
              control = list(target_basis = "cell", instrument_basis = "cell"))
  expect_lt(max(abs(predict(f, x) - y)), 1e-5)
  expect_equal(f$target_spec$type, "cell")
  expect_equal(predict(f, matrix(c(2, 2), nrow = 1)), f$coefficients[1L])
  expect_error(predict(f, matrix(1, nrow = 1)), "wrong number of columns")
})

test_that("spline learners predict finitely when a training-constant input varies later", {
  x <- cbind(seq(-1, 1, length.out = 100), 1)
  newx <- cbind(x[, 1], 0)
  y <- 1 + x[, 1]
  for (method in c("sieve_md", "landweber")) {
    expect_silent(f <- fit_cm(y, target = x, instrument = x, method = method))
    expect_true(all(is.finite(predict(f, newx))))
    expect_equal(predict(f, newx), predict(f, x))
  }
})

test_that("grouped joint-category sieve agrees with the independent dense normal equations", {
  set.seed(904)
  x <- cbind(sample(0:2, 180, TRUE), sample(0:1, 180, TRUE))
  z <- cbind(sample(0:3, 180, TRUE), x[, 2])
  d <- sample(c(0, .5, 1, -1), 180, TRUE)
  y <- rnorm(180)
  x[d == 0, ] <- NA
  f <- fit_cm(y, d, x, z, "sieve_md", list(target_basis = "cell", instrument_basis = "cell",
                                          lambda = .01, weight_ridge = .001))
  active <- d != 0
  H <- matrix(0, length(y), f$target_design_cols)
  H[active, ] <- cmbridge:::.eval_basis_spec(x[active, , drop = FALSE], f$target_spec)
  Q <- cmbridge:::.eval_basis_spec(z, f$instrument_spec)
  A <- crossprod(Q, d * H) / length(y)
  W <- cmbridge:::.safe_inverse(crossprod(Q) / length(y) + .001 * diag(ncol(Q)), ridge = .001)
  P <- diag(ncol(H)); P[1, 1] <- 0
  coef <- as.numeric(cmbridge:::.safe_solve(crossprod(A, W %*% A) + .01 * P,
                                          crossprod(A, W %*% (crossprod(Q, y) / length(y)))))
  expect_equal(f$coefficients, coef, tolerance = 1e-9)
  expect_equal(fitted(f)[active], as.numeric(H[active, , drop = FALSE] %*% coef), tolerance = 1e-9)
  expect_equal(as.matrix(f$instrument_features(z)), Q[, -1L, drop = FALSE])
  expect_equal(moment_loss(f), moment_loss(f, y, d, x, z), tolerance = 1e-10)

})

test_that("spline learners extrapolate finitely with all internal knots at a boundary", {
  for (x in list(c(1, rep(2, 8)), c(rep(1, 8), 2))) {
    x <- matrix(x, ncol = 1)
    for (method in c("sieve_md", "landweber")) {
      f <- suppressWarnings(fit_cm(1 + x[, 1], target = x, instrument = x, method = method))
      spec <- f$target_spec$specs[[1L]]
      expect_true(all(spec$knots > min(spec$Boundary.knots) & spec$knots < max(spec$Boundary.knots)))
      expect_true(all(is.finite(suppressWarnings(predict(f, matrix(c(0, 3), ncol = 1))))))
    }
  }
})

test_that("spline learners retain their prediction schema for zero measured rows", {
  x <- cbind(seq(-1, 1, length.out = 80), 1)
  empty <- x[FALSE, , drop = FALSE]
  for (method in c("sieve_md", "landweber")) {
    f <- fit_cm(1 + x[, 1], target = x, instrument = x, method = method)
    design <- cmbridge:::.eval_basis_spec(empty, f$target_spec)
    expect_equal(dim(design), c(0L, length(f$coefficients)))
    expect_identical(predict(f, empty), numeric())
    expect_equal(dim(f$instrument_features(empty)),
                 c(0L, ncol(f$instrument_features(x))))
    expect_error(predict(f, matrix(numeric(), ncol = 1)), "wrong number of columns")
  }
})
