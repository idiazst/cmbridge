# Replot saved predictions without fitting or using truth to choose weights.
plot_ensemble_truth <- function(out_dir) {
  curves <- read.csv(file.path(out_dir, "seed1_curves.csv"))
  config <- readRDS(file.path(out_dir, "configuration.rds"))
  stopifnot(all(c("scenario", "problem", "truth", "ensemble") %in% names(curves)),
    all(is.finite(curves$truth)), all(is.finite(curves$ensemble)))

  draw <- function() {
    par(mfrow = c(2, 2), mar = c(4.3, 4.5, 3.2, 1), oma = c(0, 0, 3, 0))
    for (scenario in c("cubic_truth", "rbf_truth")) {
      for (problem in c("bridge", "adjoint")) {
        cc <- curves[curves$scenario == scenario & curves$problem == problem, ]
        stopifnot(nrow(cc) > 0L)
        limits <- range(c(cc$truth, cc$ensemble))
        label <- if (scenario == "cubic_truth") "Cubic truth" else "Gaussian kernel truth"
        plot(cc$truth, cc$ensemble, pch = 19, cex = 0.55, col = "#2463A6",
          xlim = limits, ylim = limits, asp = 1,
          xlab = paste("True", problem), ylab = paste("Estimated ensemble", problem),
          main = paste(label, "-", problem))
        abline(a = 0, b = 1, lty = 2, lwd = 1.5, col = "#333333")
        legend("topleft", c("Ensemble", "y = x"),
          col = c("#2463A6", "#333333"), pch = c(19, NA),
          lty = c(NA, 2), bty = "n", cex = 0.8)
        mtext(sprintf("Grid RMSE = %.4f", sqrt(mean((cc$ensemble - cc$truth)^2))),
          side = 3, line = 0.3, cex = 0.8)
      }
    }
    mtext(sprintf("Ensemble recovery | n = %s | seed = %s",
      format(max(config$sizes), big.mark = ",", scientific = FALSE), config$seeds[1L]),
      outer = TRUE, font = 2, cex = 1.1, line = 1)
  }

  render <- function(open_device) {
    open_device()
    on.exit(dev.off())
    draw()
  }
  render(function() pdf(file.path(out_dir, "truth_vs_estimate.pdf"),
    width = 10, height = 10, bg = "white"))
  render(function() png(file.path(out_dir, "truth_vs_estimate.png"),
    width = 1600, height = 1600, res = 160))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 1L) stop("Usage: Rscript tests/validation/plot_ensemble_truth.R <ensemble-results-dir>")
  plot_ensemble_truth(args[[1L]])
}
