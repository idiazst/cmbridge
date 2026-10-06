test_that("joint cells solve a known missing-data bridge", {
  # An exact finite-frequency dataset avoids a probabilistic assertion.
  B <- V <- M <- numeric()
  count <- matrix(c(70,10,10,10,10,70,10,10,10,10,70,10,10,10,10,70),4,4)
  pi <- c(.2,.4,.6,.8)
  for(b in 1:4) for(v in 1:4) for(m in 0:1) {
    n <- round(count[b,v] * if(m) pi[v] else 1-pi[v])
    B <- c(B,rep(b,n)); V <- c(V,rep(if(m)v else NA,n)); M <- c(M,rep(m,n))
  }
  f <- fit_bridge(B,V,M,"saturated_l1",list(penalty=1e-8))
  expect_equal(predict(f,1:4),1/pi,tolerance=2e-4)
  expect_lt(f$moment_loss,1e-8)
  z <- fit_cm(rep(0,length(B)),M,V,B,"saturated_l1")
  expect_equal(predict(z,1:4),rep(0,4))
  expect_error(fit_bridge(B,V,M,"saturated_l1",list(penalty=.001,lower=1,upper=3)),
               "Post-fit prediction bounds")
})

test_that("fast cell keys preserve the original representation and row order", {
  x <- data.frame(a = rep(c(0, 1 / 6, 1, NA_real_), 100), b = rep(0:1, 200))
  expect_identical(cmbridge:::.cell_keys(as.matrix(x)), do.call(paste, c(x, sep = "|")))
})

test_that("L1 solver failures cannot silently become empty fitted models", {
  set.seed(541)
  n <- 500
  B <- sample(1:4,n,TRUE); V <- sample(1:4,n,TRUE); M <- rbinom(n,1,.5)
  V[M==0] <- NA
  expect_error(suppressWarnings(fit_bridge(B,V,M,"saturated_l1",
    list(penalty=1e-8,max_iter=1))),"failed to converge")
  expect_error(fit_bridge(B,V,M,"saturated_l1",list(penalty=-1)),"nonnegative")
})

test_that("penalty CV selects a positive scale and refits adjoint loadings within each split", {
  set.seed(701)
  n <- 600L
  B <- sample(1:4, n, TRUE); V <- sample(1:4, n, TRUE)
  M <- rbinom(n, 1, .7); V[M == 0] <- NA_real_
  folds <- sample(rep(1:2, length.out = n))
  calls <- list()
  loading <- function(train, validation) {
    expect_length(intersect(train, validation), 0L)
    calls[[length(calls) + 1L]] <<- list(train = train, valid = validation)
    list(train = rep(2, length(train)), validation = rep(2, length(validation)))
  }
  split_calls <- list()
  shared_split <- function(ids) {
    split_calls[[length(split_calls) + 1L]] <<- ids
    (seq_along(ids) - 1L) %% 2L + 1L
  }
  lib <- list(cv = list(method = "saturated_l1", control = list(
    penalty_scales = c(.001, .01, .1), penalty_ids = seq_len(n),
    penalty_folds = shared_split)))
  f <- fit_adjoint_ensemble(B, V, M, loading, lib, fold_id = folds, kernel = "cell")
  selected <- f$candidates[[1]]$penalty_cv
  expect_equal(sum(selected$selected), 1L)
  expect_true(all(is.finite(selected$loss)))
  expect_equal(f$candidates[[1]]$tuning$penalty,
               selected$scale[selected$selected] / sqrt(n))
  # Three candidate fits (two ensemble training groups plus final fit), each
  # with two penalty validation groups, in addition to ensemble loadings.
  expect_length(split_calls, 3L)
  expect_identical(split_calls[[3]], seq_len(n))
  expect_length(calls, 9L)
  expect_true(any(vapply(calls, function(x) length(x$train) < n / 2, FALSE)))
})
