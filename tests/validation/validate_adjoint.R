library(cmbridge)

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[[1L]] else file.path("results", "adjoint")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## Strong known-truth adjoint diagnostic.
##
## Let S in {-1,+1}, B = V + delta*S, and let observation depend on S:
##   P(M=1 | V,S=+1) = p_plus
##   P(M=1 | V,S=-1) = p_minus.
## Hence M depends on B beyond V, so this is not the easy MAR case.
## Among complete cases,
##   P(S=+1 | V,M=1) = p_plus/(p_plus+p_minus).
## Therefore, for any chosen true lambda,
##   phi(V) = E[lambda(B) | V,M=1]
##          = w_plus lambda(V+delta) + w_minus lambda(V-delta)
## makes
##   E[M{lambda(B)-phi(V)} | V] = 0
## hold exactly.

delta <- 0.15
p_plus <- 0.80
p_minus <- 0.40
w_plus <- p_plus / (p_plus + p_minus)
w_minus <- p_minus / (p_plus + p_minus)

generate_adjoint <- function(n, lambda_truth, seed) {
  set.seed(seed)
  V <- runif(n, -0.75, 0.75)
  S <- sample(c(-1, 1), n, replace = TRUE)
  B <- V + delta * S
  pM <- ifelse(S == 1, p_plus, p_minus)
  M <- rbinom(n, 1, pM)
  phi <- w_plus * lambda_truth(V + delta) +
    w_minus * lambda_truth(V - delta)
  list(B = B, V = V, M = M, phi = phi, S = S)
}

score_fit <- function(fit, truth, grid, method, seed) {
  est <- predict(fit, grid)
  tru <- truth(grid)
  data.frame(
    method = method,
    seed = seed,
    rmse = sqrt(mean((est - tru)^2)),
    max_abs = max(abs(est - tru)),
    moment_loss = fit$moment_loss
  )
}

summarize_results <- function(res) {
  do.call(rbind, lapply(split(res, res$method), function(x) {
    data.frame(
      method = x$method[[1L]],
      mean_rmse = mean(x$rmse),
      max_rmse = max(x$rmse),
      mean_max_abs = mean(x$max_abs),
      max_abs = max(x$max_abs),
      mean_moment_loss = mean(x$moment_loss)
    )
  }))
}

lambda_poly <- function(b) 1.2 + 0.40 * b + 0.15 * b^2 - 0.08 * b^3

rbf_truth_factory <- function(centers, bandwidth, coef) {
  force(centers); force(bandwidth); force(coef)
  function(b) as.numeric(rbf_kernel(b, centers, bandwidth) %*% coef)
}

seeds <- 1:5
grid <- seq(-0.90, 0.90, length.out = 201)
results <- list()
curves <- list()

## Sieve minimum distance: lambda is exactly cubic.
for (seed in seeds) {
  dat <- generate_adjoint(300000, lambda_poly, seed)
  fit <- fit_adjoint(
    dat$B, dat$V, dat$M, dat$phi, method = "sieve_md",
    control = list(
      target_basis = "poly", instrument_basis = "poly",
      target_degree = 3L, instrument_degree = 7L,
      lambda = 1e-10, weight_ridge = 1e-8
    )
  )
  results[[length(results) + 1L]] <-
    score_fit(fit, lambda_poly, grid, "sieve_md", seed)
  if (seed == 1L) {
    curves[["sieve_md"]] <- data.frame(
      method = "sieve_md", grid = grid,
      truth = lambda_poly(grid), estimate = predict(fit, grid)
    )
  }
}

## Landweber: same exact polynomial class.
for (seed in seeds) {
  dat <- generate_adjoint(300000, lambda_poly, seed)
  fit <- fit_adjoint(
    dat$B, dat$V, dat$M, dat$phi, method = "landweber",
    control = list(
      target_basis = "poly", instrument_basis = "poly",
      target_degree = 3L, instrument_degree = 7L,
      n_iter = 1500L, step_fraction = 0.95,
      weight_ridge = 1e-8, tol = 0
    )
  )
  results[[length(results) + 1L]] <-
    score_fit(fit, lambda_poly, grid, "landweber", seed)
  if (seed == 1L) {
    curves[["landweber"]] <- data.frame(
      method = "landweber", grid = grid,
      truth = lambda_poly(grid), estimate = predict(fit, grid)
    )
  }
}

## PMMR: lambda lies exactly in the supplied Gaussian target-RKHS/Nystrom span.
centers <- matrix(seq(-0.75, 0.75, length.out = 5), ncol = 1)
bw <- 0.65
coef <- c(0.50, 0.62, 0.55, 0.60, 0.48)
lambda_rbf <- rbf_truth_factory(centers, bw, coef)
critic_centers <- matrix(seq(-0.75, 0.75, length.out = 81), ncol = 1)

for (seed in seeds) {
  dat <- generate_adjoint(300000, lambda_rbf, seed)
  fit <- fit_adjoint(
    dat$B, dat$V, dat$M, dat$phi, method = "pmmr",
    control = list(
      scale = FALSE,
      target_centers_matrix = centers,
      instrument_centers_matrix = critic_centers,
      target_bandwidth = bw,
      instrument_bandwidth = 0.35,
      lambda = 1e-5,
      nystrom_ridge = 1e-10,
      nystrom_tol = 1e-10
    )
  )
  results[[length(results) + 1L]] <-
    score_fit(fit, lambda_rbf, grid, "pmmr", seed)
  if (seed == 1L) {
    curves[["pmmr"]] <- data.frame(
      method = "pmmr", grid = grid,
      truth = lambda_rbf(grid), estimate = predict(fit, grid)
    )
  }
}

res <- do.call(rbind, results)
summary <- summarize_results(res)
curve_df <- do.call(rbind, curves)

write.csv(res, file.path(out_dir, "validation_results.csv"), row.names = FALSE)
write.csv(summary, file.path(out_dir, "validation_summary.csv"), row.names = FALSE)
write.csv(curve_df, file.path(out_dir, "seed1_curves.csv"), row.names = FALSE)

pdf(file.path(out_dir, "truth_vs_estimate.pdf"), width = 6, height = 6)
for (nm in names(curves)) {
  cc <- curves[[nm]]
  plot(
    cc$truth, cc$estimate, pch = 19, cex = 0.45,
    xlab = "True adjoint", ylab = "Estimated adjoint", main = nm
  )
  abline(0, 1, lty = 2)
}
dev.off()

## Store a simple check that the observation mechanism really differs by S.
dat_check <- generate_adjoint(200000, lambda_poly, 991)
mechanism <- data.frame(
  group = c("S=-1", "S=+1"),
  empirical_P_M1 = c(
    mean(dat_check$M[dat_check$S == -1]),
    mean(dat_check$M[dat_check$S == 1])
  ),
  target_P_M1 = c(p_minus, p_plus)
)
write.csv(mechanism, file.path(out_dir, "observation_mechanism_check.csv"),
          row.names = FALSE)

print(res)
print(summary)
print(mechanism)

## Failure thresholds are deliberately wider than normal Monte Carlo error.
bad <- any(res$rmse >= 0.03) || any(res$max_abs >= 0.07)
if (bad) stop("Adjoint validation failed tolerance checks.", call. = FALSE)

cat("Adjoint validation passed. Outputs:", out_dir, "\n")
