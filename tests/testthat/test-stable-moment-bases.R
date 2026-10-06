test_that("nonredundant cell critics preserve the original ridge-weighted loss", {
  z <- matrix(c(rep(0, 5), rep(1, 17), rep(2, 9)), ncol = 1)
  spec <- cmbridge:::.fit_basis_spec(z, "cell")
  full <- cmbridge:::.eval_basis_spec(z, spec)
  ridge <- 1e-4
  reference <- solve(crossprod(full) / nrow(z) + ridge * diag(ncol(full)))
  reduced <- cmbridge:::.moment_basis_spec(z, spec, ridge)
  r <- cbind(seq_len(nrow(z)) / nrow(z), sin(seq_len(nrow(z))))
  m <- crossprod(full, r) / nrow(z)
  mr <- crossprod(reduced$Q, r) / nrow(z)
  expect_equal(crossprod(m, reference %*% m), crossprod(mr, reduced$W %*% mr), tolerance = 1e-10)
  expect_equal(reduced$features(z), reduced$Q)
})

test_that("large cell ridge retains the unpenalized common direction", {
  S <- crossprod(matrix(c(1,2,3,2,1,0,3,0,2),3)) / 20
  b <- c(.3,.7,.4)
  got <- cmbridge:::.centered_ridge_solve(S,b,1e14)
  limit <- sum(b)/sum(S)
  expect_equal(got,rep(limit,3),tolerance=1e-10)
  expect_equal(sum(S %*% got-b),0,tolerance=1e-10)
  # At a moderate penalty compare the independent dense equations.
  center <- diag(3)-matrix(1/3,3,3)
  expect_equal(cmbridge:::.centered_ridge_solve(S,b,.1),as.numeric(solve(S+.1*center,b)),tolerance=1e-10)
})

test_that("joint cells retain fixed main effects for unseen combinations", {
  x <- as.matrix(expand.grid(a=0:1,b=0:1))
  spec <- cmbridge:::.fit_basis_spec(x[-4,,drop=FALSE],"cell_linear")
  H <- cmbridge:::.eval_basis_spec(x,spec)
  theta <- c(.25,-.4,.6,rep(0,length(spec$levels)))
  expect_equal(as.numeric(H%*%theta),.25-.4*x[,1]+.6*x[,2])
  expect_identical(ncol(H),1L+ncol(x)+length(spec$levels))
  target <- x[rep(1:4,each=50),]
  response <- 1+exp(-(.25-.4*target[,1]+.6*target[,2]))
  f <- fit_cm(response,target=target,instrument=target,method="sieve_md",
    control=list(target_basis="cell_linear",instrument_basis="cell",link="inverse_logit",lambda=1e-8))
  expect_lt(max(abs(predict(f,x)-(1+exp(-(.25-.4*x[,1]+.6*x[,2]))))),1e-3)
  expect_equal(moment_loss(f),moment_loss(f,response,rep(1,length(response)),target,target),tolerance=1e-10)
})
