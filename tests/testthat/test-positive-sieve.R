test_that("inverse-logit saturated sieve contains positive interaction functions", {
  grid <- expand.grid(a = 0:1, b = 0:1)
  x <- as.matrix(grid[rep(1:4, each = 100), ])
  y <- 1.2 + 2 * x[, 1] * x[, 2] + x[, 2]
  f <- fit_cm(y, target = x, instrument = x, method = "sieve_md", control = list(
    target_basis = "cell", instrument_basis = "cell", link = "inverse_logit"))
  expect_lt(max(abs(predict(f, x) - y)), 1e-5)
  expect_true(all(predict(f, x) >= 1))
  expect_equal(predict(f, x), 1 + exp(-f$link_coefficients[match(
    cmbridge:::.cell_keys(x), f$target_spec$levels)]), tolerance = 1e-12)
  expect_gte(predict(f, matrix(c(2, 2), 1)), 1)
  expect_identical(predict(f, x[FALSE, , drop = FALSE]), numeric())
  expect_true(f$solver$converged)
  expect_equal(moment_loss(f), moment_loss(f, y, diagonal = rep(1, length(y)), target = x, instrument = x), tolerance = 1e-10)
  expect_equal(predict(unserialize(serialize(f, NULL)), x), predict(f, x))
})

test_that("unrestricted inverse-expit coefficients optimize the moment loss", {
  S <- matrix(c(2, 1, 1, 2), 2)
  b <- c(0, 5)
  constant <- as.numeric(crossprod(b, solve(S, b)))
  problem <- cmbridge:::.sieve_link_objective(S, b, constant, "inverse_logit")
  par <- c(.3, -.7)
  eps <- 1e-5
  numeric_gradient <- vapply(1:2, function(j) {
    step <- numeric(2); step[j] <- eps
    (problem$fn(par + step) - problem$fn(par - step)) / (2 * eps)
  }, numeric(1))
  expect_equal(problem$gr(par), numeric_gradient, tolerance = 1e-6)
  got <- cmbridge:::.sieve_link_solve(S, b, constant, "inverse_logit", 1e-10, 10000)
  expect_true(all(is.finite(got$parameters)))
  expect_true(got$converged)
  expect_equal(got$values, 1 + exp(-got$parameters))
  # Independent analytic infimum: beta=(1,2), with first coefficient tending
  # to +Inf. The finite unconstrained coefficients approach that solution.
  expect_equal(got$values, c(1, 2), tolerance = 1e-4)
  loss <- function(v) as.numeric(crossprod(v, S %*% v) - 2 * crossprod(b, v))
  clipped <- pmax(1, solve(S, b))
  expect_lt(loss(got$values), loss(clipped) - 1)
  expect_true(is.null(got$active_constraints))
})

test_that("finite unrestricted coefficients approach beta=1 and adjoints stay unrestricted", {
  x <- matrix(rep(0:1, each = 30), ncol = 1)
  controls <- list(target_basis = "cell", instrument_basis = "cell", link = "inverse_logit", lambda = .1)
  f <- fit_cm(rep(-1, 60), target = x, instrument = x, method = "sieve_md", control = controls)
  expect_equal(predict(f, x), rep(1, 60), tolerance = 1e-4)
  expect_true(all(is.finite(f$link_coefficients)))
  expect_equal(predict(f, matrix(2, 1)), 1, tolerance = 1e-4)
  controls$link <- "identity"
  adj <- fit_adjoint(x, x, rep(1, 60), rep(-1, 60), "sieve_md", controls)
  expect_true(all(predict(adj, x) < 0))
  controls$link <- "log"
  positive <- fit_cm(rep(-1, 60), target = x, instrument = x, method = "sieve_md", control = controls)
  expect_true(all(predict(positive, x) > 0))
  expect_equal(positive$residual, rep(-1, 60) - fitted(positive))
  expect_error(fit_cm(rep(1, 60), target = x, instrument = x, method = "sieve_md",
    control = list(link = "inverse_logit")), "require.*cell")
})

test_that("loss evaluation remains stable near an inverse-expit boundary", {
  # Regression case from a small nested training fit. Subtracting the original
  # large quadratic terms at every evaluation obscured changes near zero loss.
  S <- matrix(c(.0250000066385266,.0249999972457489,-2.16269469430132e-9,-2.08472878622389e-9,
    .0249999972457489,.125000004277624,-2.09673732579113e-9,-1.9644895189369e-9,
    -2.16269469430132e-9,-2.09673732579113e-9,.0233333402890839,.00999999733084854,
    -2.08472878622389e-9,-1.9644895189369e-9,.00999999733084854,.0500000054138374),4)
  b <- c(.0999999986235052,.199999997839332,.0999999990146607,.099999998388812)
  got <- cmbridge:::.sieve_link_solve(S,b,.99999998800084,"inverse_logit",1e-12,10000)
  expect_true(all(is.finite(got$parameters)))
  expect_true(got$converged)
  expect_equal(got$values,as.numeric(solve(S,b)),tolerance=1e-4)
  expect_lte(got$function_gradient,got$tolerance)
})
test_that("penalty CV flat-link coefficients are refined without relaxing function accuracy", {
  case <- readRDS(test_path("fixtures", "sieve-cv-flat-link.rds"))
  fit <- cmbridge:::.sieve_link_solve(case$S, case$b, case$constant, case$link, 1e-10, 10000L)
  expect_true(fit$converged)
  expect_true(all(is.finite(fit$parameters)))
  expect_true(all(fit$values >= 1))
  expect_lte(max(fit$function_gradient, fit$coefficient_gradient), fit$tolerance)
  # Independently evaluate the stationary function values with coordinate 2
  # at the limiting value; the production optimizer has no coefficient bounds.
  free <- c(1, 3, 4)
  reference <- rep(1, 4)
  reference[free] <- solve(case$S[free, free], case$b[free] - case$S[free, 2])
  expect_equal(fit$values, reference, tolerance = 1e-4)
})

test_that("link-step scaling is updated after a large unrestricted coefficient move", {
  case <- readRDS(test_path("fixtures", "sieve-cv-refinement.rds"))
  fit <- cmbridge:::.sieve_link_solve(case$S, case$b, case$constant, case$link, 1e-10, 10000L)
  expect_true(fit$converged)
  expect_true(all(is.finite(fit$parameters)))
  expect_lte(max(fit$function_gradient, fit$coefficient_gradient), fit$tolerance)
  free <- c(1, 3)
  limiting <- c(2, 4)
  reference <- rep(1, 4)
  reference[free] <- solve(case$S[free, free], case$b[free] - rowSums(case$S[free, limiting]))
  expect_equal(fit$values, reference, tolerance = 1e-4)
})

test_that("small nested CV link fits retain function accuracy near stationarity", {
  case <- readRDS(test_path("fixtures", "sieve-cv-nested-precision.rds"))
  fit <- cmbridge:::.sieve_link_solve(case$S, case$b, case$constant, case$link, 1e-10, 10000L)
  expect_true(fit$converged)
  expect_true(all(is.finite(fit$parameters)))
  expect_lte(max(fit$function_gradient, fit$coefficient_gradient), fit$tolerance)
  expect_true(all(fit$values >= 1))
})

test_that("a revived inverse-link coefficient is refined after its restart", {
  case <- readRDS(test_path("fixtures", "sieve-cv-restart-refinement.rds"))
  fit <- cmbridge:::.sieve_link_solve(case$S, case$b, case$constant, case$link, 1e-10, 10000L)
  expect_true(fit$converged)
  expect_true(all(is.finite(fit$parameters)))
  expect_lte(max(fit$function_gradient, fit$coefficient_gradient), fit$tolerance)
  expect_true(any(grepl("restart", fit$trace$method)))
  expect_lt(tail(fit$trace$coefficient_gradient, 1), 1e-5)
})
