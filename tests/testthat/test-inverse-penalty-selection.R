test_that("ridge GCV uses the inverse smoother trace and native residual geometry", {
  set.seed(610)
  z <- cbind(rnorm(140), rnorm(140))
  x <- cbind(z[, 1] + rnorm(140), .1 * z[, 2] + rnorm(140))
  y <- 1 + x[, 1] - 2 * x[, 2] + rnorm(140)
  for (method in c("sieve_md", "pmmr")) {
    control <- if (method == "sieve_md") list(target_basis = "poly", target_degree = 1,
      instrument_basis = "quadratic", lambda = .01) else list(lambda = .01, target_centers = 8, instrument_centers = 12)
    fit <- fit_cm(y, target = x, instrument = z, method = method, control = control)
    g <- cmbridge:::.inverse_geometry(fit, y, rep(1, 140), x, z)
    expected <- sum(diag(g$J %*% solve(crossprod(g$J) + .01 * g$P, t(g$J))))
    expect_equal(cmbridge:::.inverse_ridge_df(g, fit, rep(1, 140)), expected, tolerance = 1e-7)
    expect_equal(g$native_loss, fit$moment_loss, tolerance = 1e-8)
  }
})

test_that("nonlinear sieve sensitivity matches perturbations of the native optimum", {
  set.seed(611)
  z <- matrix(rnorm(180), ncol = 1)
  x <- z + matrix(rnorm(180), ncol = 1)
  d <- rbinom(180, 1, .65); y <- rep(1, 180)
  x[d == 0, ] <- NA
  ctrl <- list(target_basis = "poly", target_degree = 1,
    instrument_basis = "quadratic", link = "inverse_logit", penalty_on = "link",
    lambda = .005, solver_tolerance = 1e-14)
  fit <- fit_cm(y, d, x, z, "sieve_md", ctrl)
  g <- cmbridge:::.inverse_geometry(fit, y, d, x, z)
  lift <- crossprod(g$L, solve(tcrossprod(g$L)))
  epsilon <- 1e-4
  calculated <- 0
  for (j in seq_len(g$rank)) {
    fits <- lapply(c(-1, 1), function(sign) fit_cm(y + sign * epsilon * lift[, j], d, x, z, "sieve_md", ctrl))
    delta <- g$L %*% (d * (ifelse(is.na(fits[[2]]$fitted), 0, fits[[2]]$fitted) -
      ifelse(is.na(fits[[1]]$fitted), 0, fits[[1]]$fitted))) / (2 * epsilon)
    calculated <- calculated + delta[j]
  }
  expect_equal(cmbridge:::.inverse_ridge_df(g, fit, d), as.numeric(calculated), tolerance = 1e-3)
})

test_that("Landweber trace includes its data-dependent initial intercept", {
  set.seed(612)
  z <- matrix(rnorm(120), ncol = 1); x <- z + matrix(rnorm(120), ncol = 1)
  y <- 1 + x[, 1] + rnorm(120); d <- rep(1, 120)
  ctrl <- list(target_basis = "poly", target_degree = 1, instrument_basis = "quadratic", n_iter = 12)
  fit <- fit_cm(y, d, x, z, "landweber", ctrl)
  g <- cmbridge:::.inverse_geometry(fit, y, d, x, z)
  options <- cmbridge:::.inverse_selection_control(list())
  actual <- cmbridge:::.inverse_landweber_df(g, fit, y, d, x, z, ctrl, 1 / 12, options)
  lift <- crossprod(g$L, solve(tcrossprod(g$L))); calculated <- 0
  for (j in seq_len(g$rank)) {
    fits <- lapply(c(-1, 1), function(sign) fit_cm(y + sign * 1e-5 * lift[, j], d, x, z, "landweber", ctrl))
    calculated <- calculated + as.numeric(g$L %*% (fits[[2]]$fitted - fits[[1]]$fitted))[j] / 2e-5
  }
  expect_equal(actual$df, calculated, tolerance = 1e-6)
})

test_that("L-curve selects positive corner curvature in a known inverse path", {
  A <- rbind(diag(c(1, .1, .01)), matrix(0, 3, 3))
  b <- c(1, .12, .04, .005, -.005, .005)
  grid <- 10^seq(-7, 1, length.out = 81)
  records <- lapply(grid, function(lambda) {
    theta <- solve(crossprod(A) + lambda * diag(3), crossprod(A, b))
    list(residual_norm = sqrt(sum((A %*% theta - b)^2)), seminorm = sqrt(sum(theta^2)), values = theta, gcv = 0)
  })
  scores <- cmbridge:::.inverse_path_scores(grid, records, "lcurve", cmbridge:::.inverse_selection_control(list()))
  choice <- which.min(scores$score)
  expect_true(choice > 1 && choice < length(grid))
  expect_gt(scores$curvature[choice], 0)
  sigma2 <- c(1, .1, .01)^2; signal2 <- b[1:3]^2
  analytic <- vapply(grid, function(lambda) {
    r <- sum(signal2 * (lambda / (sigma2 + lambda))^2) + sum(b[4:6]^2)
    rp <- sum(2 * signal2 * lambda * sigma2 / (sigma2 + lambda)^3)
    rpp <- sum(2 * signal2 * sigma2 * (sigma2 - 2 * lambda) / (sigma2 + lambda)^4)
    v <- sum(sigma2 * signal2 / (sigma2 + lambda)^2)
    vp <- -2 * sum(sigma2 * signal2 / (sigma2 + lambda)^3)
    vpp <- 6 * sum(sigma2 * signal2 / (sigma2 + lambda)^4)
    dx <- .5 * lambda * rp / r; dy <- .5 * lambda * vp / v
    ddx <- .5 * (lambda * rp / r + lambda^2 * rpp / r - (lambda * rp / r)^2)
    ddy <- .5 * (lambda * vp / v + lambda^2 * vpp / v - (lambda * vp / v)^2)
    (dy * ddx - dx * ddy) / (dx^2 + dy^2)^1.5
  }, 0)
  expect_lte(abs(log(grid[choice] / grid[which.max(analytic)])), .25)
  expect_equal(scores$curvature[choice], analytic[choice], tolerance = .1)
})

test_that("quasi-optimality cannot accept a constant path", {
  grid <- 10^seq(-6, 0, length.out = 13)
  records <- lapply(grid, function(lambda) list(residual_norm = 1, seminorm = 1, values = rep(2, 10), gcv = 0))
  scores <- cmbridge:::.inverse_path_scores(grid, records, "quasiopt", cmbridge:::.inverse_selection_control(list()))
  expect_true(all(!is.finite(scores$score)))
})

test_that("selected fits or explicit failures preserve positive paths and RNG", {
  set.seed(613)
  z <- matrix(rnorm(90), ncol = 1); x <- z + matrix(rnorm(90), ncol = 1)
  y <- 1 + x[, 1] + rnorm(90)
  for (method in c("sieve_md", "landweber", "pmmr")) for (rule in c("gcv", "lcurve", "quasiopt")) {
    ctrl <- if (method == "pmmr") list(target_centers = 6, instrument_centers = 10) else
      list(target_basis = "poly", target_degree = 1, instrument_basis = "quadratic")
    state <- .Random.seed
    answer <- tryCatch(select_cm_penalty(y, target = x, instrument = z, method = method, control = ctrl,
      rule = rule, selection_control = list(scales = 10^seq(-6, 0, length.out = 13),
        iteration_grid = c(2, 4, 8, 16, 32, 64), max_extensions = 1)), error = identity)
    expect_identical(.Random.seed, state)
    if (inherits(answer, "error")) {
      expect_s3_class(answer, "cmbridge_penalty_selection_error")
      expect_true(all(answer$penalty_path$scale > 0))
      expect_gt(nrow(answer$boundary_history), 0)
    } else {
      expect_identical(answer$penalty_selection, rule)
      at <- which(answer$penalty_path$selected)
      expect_true(at > 1 && at < nrow(answer$penalty_path))
      expect_length(at, 1)
      expect_true(all(is.finite(predict(answer, x))))
      expect_identical(answer$penalty_training_n, length(y))
    }
  }
})


test_that("pure-cell bridge geometry respects function versus link regularization", {
  set.seed(615)
  x <- matrix(rbinom(240, 1, .5), ncol = 1)
  z <- matrix((x[, 1] + rbinom(240, 1, .45)) %% 3L, ncol = 1)
  d <- rbinom(240, 1, .55 + .15 * x[, 1]); y <- rep(1, 240)
  x[d == 0, ] <- NA
  for (penalty in c("function", "link")) {
    fit <- fit_cm(y, d, x, z, "sieve_md", list(target_basis = "cell", instrument_basis = "cell",
      link = "inverse_logit", penalty_on = penalty, lambda = .01))
    geometry <- cmbridge:::.inverse_geometry(fit, y, d, x, z)
    expect_equal(geometry$native_loss, fit$moment_loss, tolerance = 1e-8)
    expect_equal(geometry$values, fit$fitted[d != 0], tolerance = 1e-8)
    expect_true(is.finite(cmbridge:::.inverse_ridge_df(geometry, fit, d)))
  }
})
