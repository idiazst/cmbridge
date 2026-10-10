args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
root <- args[1L]
for (kind in c("bridge", "adjoint")) {
  folder <- file.path(root, kind)
  curves <- read.csv(file.path(folder, "seed1_curves.csv"))
  draw <- function() {
    par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))
    for (method in c("sieve_md", "landweber", "pmmr")) {
      d <- curves[curves$method == method, ]
      limits <- range(d$truth, d$estimate)
      plot(d$truth, d$estimate, xlim = limits, ylim = limits, asp = 1,
        xlab = "True value", ylab = "Estimated value", pch = 19, cex = .5,
        col = "#2166ac", main = paste(kind, method))
      abline(0, 1, lty = 2, col = "#555555")
    }
  }
  png(file.path(folder, "truth_vs_estimate.png"), width = 1800, height = 600, res = 150)
  draw(); dev.off()
}
