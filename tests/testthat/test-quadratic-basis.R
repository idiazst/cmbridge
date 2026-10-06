test_that("the complete quadratic basis is fixed and includes every pairwise interaction", {
  x <- as.matrix(expand.grid(a = 0:3, b = 0:1, c = 0:1))
  spec <- cmbridge:::.fit_basis_spec(x, "quadratic")
  expect_identical(spec, cmbridge:::.fit_basis_spec(x[1, , drop = FALSE], "quadratic"))
  design <- cmbridge:::.eval_basis_spec(x, spec)
  expected <- cbind(1, x[,1], x[,1]^2, x[,2], x[,2]^2, x[,3], x[,3]^2,
    x[,1]*x[,2], x[,1]*x[,3], x[,2]*x[,3])
  expect_equal(unname(design), unname(expected))
  expect_equal(dim(cmbridge:::.eval_basis_spec(x[FALSE, , drop = FALSE], spec)), c(0L, 10L))
  # Includes interactions at combinations absent from the basis construction.
  eta <- .2 + .3*x[,1] - .4*x[,2] + .7*x[,2]*x[,3]
  theta <- c(.2,.3,0,-.4,0,0,0,0,0,.7)
  expect_equal(as.numeric(design %*% theta), eta, tolerance = 1e-12)
})
