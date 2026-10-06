test_that("grouped nonlinear moments and gradients equal the full observation calculation", {
  x <- matrix(rep(0:2, each=20), ncol=1)
  z <- cbind(1, as.numeric(x[,1]==1), as.numeric(x[,1]==2))
  d <- rep(c(0,1,1,1),15); y <- rep(c(.2,1,2),20)
  active <- d!=0; x[!active,] <- NA
  H <- cbind(1,x[active,1]); W <- diag(c(1,.7,1.3))
  design <- cmbridge:::.linked_moment_design(y,d,x,H,z,W)
  for (link in c("inverse_logit","log")) {
    problem <- cmbridge:::.linked_moment_problem(design,link,penalty=.03)
    theta <- c(.1,-.2)
    value <- cmbridge:::.cm_link(as.numeric(H%*%theta),link)
    prediction <- numeric(length(y));prediction[active]<-value
    moment <- as.numeric(crossprod(z,y-d*prediction)/length(y))
    expect_equal(problem$fn(theta),as.numeric(crossprod(moment,W%*%moment))+.03*sum(theta^2),tolerance=1e-12)
    epsilon <- 1e-5
    reference <- vapply(1:2,function(j) {
      change <- numeric(2);change[j]<-epsilon
      (problem$fn(theta+change)-problem$fn(theta-change))/(2*epsilon)
    },numeric(1))
    expect_equal(problem$gr(theta),reference,tolerance=1e-6)
  }
})

test_that("all three bridge learners use inverse-expit during estimation", {
  grid <- expand.grid(a=0:1,b=0:1)
  x <- as.matrix(grid[rep(1:4,each=300),])
  M <- rep(c(rep(1,150),rep(0,150),rep(1,120),rep(0,180),
             rep(1,100),rep(0,200),rep(1,200),rep(0,100)),1)
  observed <- x;observed[M==0,]<-NA
  truth <- c(2,2.5,3,1.5)
  controls <- list(
    sieve_md=list(target_basis="cell",instrument_basis="cell",link="inverse_logit"),
    landweber=list(target_basis="cell",instrument_basis="cell",link="inverse_logit",n_iter=5000),
    pmmr=list(link="inverse_logit",lambda=1e-8,solver_tolerance=1e-12,
      target_centers_matrix=as.matrix(grid),instrument_centers_matrix=as.matrix(grid),
      target_bandwidth=1,instrument_bandwidth=1))
  for (method in names(controls)) {
    f <- fit_bridge(x,observed,M,method,controls[[method]])
    prediction <- predict(f,as.matrix(grid))
    expect_true(all(is.finite(f$link_coefficients)))
    expect_true(all(prediction>=1))
    expect_lt(max(abs(prediction-truth)),.02)
    expect_equal(predict(f,observed[M==1,,drop=FALSE]),fitted(f)[M==1],tolerance=1e-10)
    expect_equal(moment_loss(f),moment_loss(f,rep(1,length(M)),M,observed,x),tolerance=1e-10)
    expect_equal(predict(unserialize(serialize(f,NULL)),as.matrix(grid)),prediction)
    expect_identical(predict(f,x[FALSE,,drop=FALSE]),numeric())
    if(method=="landweber") {
      expect_lte(f$solver$final_loss,f$solver$initial_loss)
      expect_lte(f$solver$final_step,f$solver$maximum_step)
    }
    if(method=="pmmr") expect_lte(f$solver$penalized_loss,f$solver$initial_loss)
  }
})

test_that("adjoints retain their unrestricted default parameterizations", {
  x <- matrix(rep(0:1,each=50),ncol=1)
  for(method in c("sieve_md","landweber","pmmr")) {
    fit <- fit_adjoint(x,x,rep(1,100),rep(-1,100),method)
    expect_identical(fit$tuning$link,"identity")
    expect_true(all(predict(fit,x)<0))
  }
})

test_that("linked Landweber uses a constant basis continuation outside training support", {
  x <- matrix(c(rep(1, 5), rep(2, 25)), ncol = 1)
  M <- rep(c(1, 0, 1), 10)
  x[M == 0, ] <- NA_real_
  f <- suppressWarnings(fit_bridge(rep(0:1, 15), x, M, "landweber",
    list(link = "inverse_logit")))
  newx <- matrix(c(-100, 1, 2, 100), ncol = 1)
  prediction <- suppressWarnings(predict(f, newx))
  expect_true(all(is.finite(prediction)))
  expect_equal(prediction[1], prediction[2])
  expect_equal(prediction[3], prediction[4])
  expect_identical(f$target_spec$extrapolation, "constant")
  identity <- suppressWarnings(fit_adjoint(rep(0:1, 15), x, M, rep(-1, 30), "landweber"))
  expect_identical(identity$tuning$target_extrapolation, "polynomial")
})

test_that("linked Landweber remains finite on the failing nested numerical-dose split", {
  x <- readRDS(test_path("fixtures", "linked-landweber-small-split.rds"))
  for (fold in unique(x$fold_id)) {
    train <- which(x$fold_id != fold); validation <- which(x$fold_id == fold & x$M == 1)
    fit <- suppressWarnings(fit_bridge(x$B[train, , drop = FALSE], x$V[train, , drop = FALSE],
      x$M[train], "landweber", list(link = "inverse_logit")))
    prediction <- predict(fit, x$V[validation, , drop = FALSE])
    expect_true(all(is.finite(prediction)))
    expect_true(all(prediction >= 1))
    expect_lt(max(prediction), 10)
    expect_lte(fit$solver$final_step, fit$solver$maximum_step)
  }
})
