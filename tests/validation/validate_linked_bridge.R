library(cmbridge)

# Supplemental finite-support recovery test for the new parameterizations.
# The existing continuous polynomial/kernel recovery tests remain unchanged.
args <- commandArgs(TRUE)
out <- if (length(args)) args[1L] else "results/linked-bridge"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
grid <- as.matrix(expand.grid(a = 0:1, b = 0:1))
beta <- c(2, 2.5, 3, 1.5)
lambda <- c(-.5, .25, 1.5, -1)
controls <- list(
  sieve_md = list(target_basis = "cell", instrument_basis = "cell", link = "inverse_logit"),
  landweber = list(target_basis = "cell", instrument_basis = "cell", link = "inverse_logit", n_iter = 5000L),
  pmmr = list(link = "inverse_logit", lambda = 1e-8, solver_tolerance = 1e-12,
    scale = FALSE, target_centers_matrix = grid, instrument_centers_matrix = grid,
    target_bandwidth = 1, instrument_bandwidth = 1))
results <- curves <- list()
for (seed in 1:5) {
  set.seed(seed)
  n <- 300000L
  cell <- sample.int(4L, n, replace = TRUE)
  B <- grid[cell, , drop = FALSE]
  M <- rbinom(n, 1, 1 / beta[cell])
  V <- B; V[M == 0, ] <- NA_real_
  for (method in names(controls)) {
    bridge <- fit_bridge(B, V, M, method, controls[[method]])
    adjoint_control <- controls[[method]]; adjoint_control$link <- "identity"
    adjoint <- fit_adjoint(B, V, M, lambda[cell], method, adjoint_control)
    for (kind in c("bridge", "adjoint")) {
      fit <- if (kind == "bridge") bridge else adjoint
      truth <- if (kind == "bridge") beta else lambda
      estimate <- predict(fit, grid)
      stopifnot(all(is.finite(estimate)),
        identical(fit$tuning$link, if (kind == "bridge") "inverse_logit" else "identity"))
      if (kind == "bridge") stopifnot(all(estimate >= 1), all(is.finite(fit$link_coefficients)))
      else stopifnot(any(estimate < 0))
      # Recovering the original conditional moments confirms that the link
      # is applied during fitting, rather than appended to an identity fit.
      stopifnot(abs(moment_loss(fit) - moment_loss(fit, fit$response, fit$diagonal,
        if (kind == "bridge") V else B[M == 1, , drop = FALSE],
        if (kind == "bridge") B else V[M == 1, , drop = FALSE])) < 1e-10)
      results[[length(results) + 1L]] <- data.frame(kind, method, seed, n,
        link = fit$tuning$link, rmse = sqrt(mean((estimate - truth)^2)),
        max_abs = max(abs(estimate - truth)), moment_loss = moment_loss(fit))
      if (seed == 1L) curves[[length(curves) + 1L]] <- data.frame(kind, method,
        cell = seq_len(4), truth, estimate)
    }
  }
}
results <- do.call(rbind, results); curves <- do.call(rbind, curves)
write.csv(results, file.path(out, "validation_results.csv"), row.names = FALSE)
summary <- do.call(rbind, lapply(split(results, list(results$kind, results$method)), function(x)
  data.frame(kind = x$kind[1], method = x$method[1], mean_rmse = mean(x$rmse),
    max_rmse = max(x$rmse), max_abs = max(x$max_abs))))
write.csv(summary, file.path(out, "validation_summary.csv"), row.names = FALSE)
write.csv(curves, file.path(out, "seed1_curves.csv"), row.names = FALSE)
draw <- function() {
  par(mfrow = c(2, 3), mar = c(4.5, 4.5, 3, 1))
  for (kind in c("bridge", "adjoint")) for (method in names(controls)) {
    z <- curves[curves$kind == kind & curves$method == method, ]
    limits <- range(z$truth, z$estimate)
    plot(z$truth, z$estimate, xlim = limits, ylim = limits, asp = 1,
      xlab = "True value", ylab = "Estimated value", pch = 19, col = "#2166ac",
      main = paste(kind, method))
    abline(0, 1, lty = 2, col = "#555555")
  }
}
pdf(file.path(out, "truth_vs_estimate.pdf"), width = 11, height = 7); draw(); dev.off()
png(file.path(out, "truth_vs_estimate.png"), width = 1650, height = 1050, res = 150); draw(); dev.off()
print(summary, row.names = FALSE)
stopifnot(all(results$rmse < .03), all(results$max_abs < .07))
cat("Inverse-expit bridge and unrestricted signed-adjoint checks passed.\n")
