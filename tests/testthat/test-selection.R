test_that("cell U-statistic matches explicit ordered pairs and retains treatment", {
  z <- cbind(H = 0, A = c(0, 0, 1, 1, 1, 2))
  r <- cbind(c(1, -1, 2, 3, 4, 99), c(3, 1, 0, 0, 1, -99))
  direct <- matrix(0, 2, 2)
  for (group in split(seq_len(nrow(z)), z[, "A"])) if (length(group) > 1) {
    for (i in group) for (j in setdiff(group, i))
      direct <- direct + tcrossprod(r[i, ], r[j, ]) / (nrow(z) * (length(group) - 1))
  }
  got <- cell_moment_gram(r, z)
  expect_equal(as.vector(got), as.vector(direct), tolerance = 1e-12)
  expect_equal(attr(got, "cell_scoring")$singleton_cells, 1)
  expect_equal(attr(got, "cell_scoring")$paired_rows, 5)
  # Removing A would hide this conditional-moment violation.
  residual <- c(1, 1, -1, -1)
  expect_equal(as.numeric(cell_moment_gram(residual, cbind(H = 0, A = c(0, 0, 1, 1)))), 1)
  expect_lt(as.numeric(cell_moment_gram(residual, matrix(0, 4, 1))), 0)
  empty_pairs <- cell_moment_gram(cbind(1:3, c(4, 5, 6)), matrix(1:3))
  expect_equal(as.vector(empty_pairs), rep(0, 4))
  expect_equal(attr(empty_pairs, "cell_scoring")$paired_rows, 0)
  expect_equal(attr(empty_pairs, "cell_scoring")$singleton_cells, 3)
})

test_that("removing self-products removes candidate-dependent noise in expectation", {
  # Enumerate all independent Rademacher errors: an exact expectation check.
  errors <- as.matrix(expand.grid(rep(list(c(-1, 1)), 4)))
  z <- matrix(c(0, 0, 1, 1), ncol = 1)
  scores <- lapply(seq_len(nrow(errors)), function(i)
    cell_moment_gram(cbind(5 * errors[i, ], .5 + errors[i, ]), z))
  expected <- Reduce(`+`, scores) / nrow(errors)
  expect_equal(as.vector(expected), as.vector(matrix(c(0, 0, 0, .25), 2)), tolerance = 1e-12)
})

test_that("PSD projection preserves the raw score and stabilizes the simplex objective", {
  raw <- matrix(c(-2, .5, .5, 1), 2)
  p <- cmbridge:::.project_gram_psd(raw)
  expect_equal(p$removed_negative_eigenvalues, 1)
  expect_gte(min(eigen(p$gram, symmetric = TRUE)$values), -1e-12)
  fit <- cmbridge:::.simplex_qp(p$gram)
  grid <- seq(0, 1, length.out = 10001)
  values <- vapply(grid, function(x) as.numeric(crossprod(c(x, 1-x), p$gram %*% c(x, 1-x))), 0)
  expect_lte(as.numeric(crossprod(fit$weights, p$gram %*% fit$weights)), min(values) + 1e-7)
})

test_that("penalty selection records boundary failures and extends both directions", {
  for (target in c(-8, 2)) {
    selected <- select_penalty_grid(c(.001, .01, .1), function(scales) (log10(scales) - target)^2)
    expect_equal(log10(selected$scale), target, tolerance = 1e-12)
    index <- which(selected$penalty_cv$selected)
    expect_true(index > 1 && index < nrow(selected$penalty_cv))
    expect_true(all(selected$boundary_history$status == "boundary_failure"))
    expect_true(all(selected$penalty_cv$scale > 0))
  }
  failure <- tryCatch(select_penalty_grid(c(.001, .01, .1), identity, max_extensions = 2), error = identity)
  expect_s3_class(failure, "cmbridge_penalty_boundary_error")
  expect_equal(nrow(failure$boundary_history), 3)
  expect_false(any(failure$penalty_cv$selected))
  tied <- select_penalty_grid(c(.001, .01, .1), function(x) rep(1, length(x)))
  expect_equal(tied$scale, .01)
})
