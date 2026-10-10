library(cmbridge)

# Only the named correct candidate contains the truth in each scenario.
# This tests selection of candidate specifications, not universal superiority
# of one estimation method. Truth is used for evaluation, never for stacking.
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[[1L]] else file.path("results", "ensemble")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
sizes <- as.integer(strsplit(Sys.getenv("CMBRIDGE_ENSEMBLE_SIZES", "2000,20000,200000"), ",")[[1L]])
seeds <- as.integer(strsplit(Sys.getenv("CMBRIDGE_ENSEMBLE_SEEDS", "1,2,3,4,5"), ",")[[1L]])
stopifnot(all(is.finite(sizes)), all(sizes > 0), all(is.finite(seeds)))
n_folds <- 3L

poly_control <- function(degree) list(target_basis = "poly", instrument_basis = "poly",
  target_degree = degree, instrument_degree = 7L, lambda = 1e-10, weight_ridge = 1e-8)
landweber_control <- function(degree) {
  ctrl <- poly_control(degree)
  ctrl$lambda <- NULL
  c(ctrl, list(n_iter = 1500L, step_fraction = 0.95, tol = 0))
}
kernel_control <- function(centers, bw) list(scale = FALSE,
  target_centers_matrix = matrix(centers, ncol = 1L), target_bandwidth = bw,
  instrument_centers_matrix = matrix(seq(-1, 1, length.out = 41), ncol = 1L),
  instrument_bandwidth = 0.30, lambda = 1e-8,
  nystrom_ridge = 1e-10, nystrom_tol = 1e-10)

rbf_centers <- seq(-0.8, 0.8, length.out = 5)
rbf_bandwidth <- 0.22
rbf_coef <- c(3.5, 1.4, 5.2, 1.4, 3.5)
rbf_truth <- function(x) as.numeric(rbf_kernel(x, rbf_centers, rbf_bandwidth) %*% rbf_coef)
scenarios <- list(
  cubic_truth = list(
    beta = function(v) 3 + 0.6 * v + 0.6 * v^2 + 1.1 * v^3,
    lambda = function(b) 1.2 + 0.4 * b + 0.6 * b^2 + 1.1 * b^3,
    correct = "sieve_cubic",
    library = list(
      sieve_cubic = list(method = "sieve_md", control = poly_control(3L)),
      landweber_linear = list(method = "landweber", control = landweber_control(1L)),
      pmmr_one_center = list(method = "pmmr", control = kernel_control(0, 0.65))
    )
  ),
  rbf_truth = list(
    beta = rbf_truth, lambda = rbf_truth, correct = "pmmr_five_centers",
    library = list(
      sieve_linear = list(method = "sieve_md", control = poly_control(1L)),
      landweber_quadratic = list(method = "landweber", control = landweber_control(2L)),
      pmmr_five_centers = list(method = "pmmr", control = kernel_control(rbf_centers, rbf_bandwidth))
    )
  )
)

generate_bridge <- function(n, truth, seed) {
  set.seed(seed)
  B <- runif(n, -1, 1)
  V <- 0.97 * B + 0.03 * runif(n, -1, 1)
  beta <- truth(V)
  stopifnot(min(beta) > 1)
  M <- rbinom(n, 1, 1 / beta)
  V[M == 0] <- NA_real_
  list(B = B, V = V, M = M)
}
generate_adjoint <- function(n, truth, seed) {
  set.seed(seed)
  V <- runif(n, -0.75, 0.75)
  S <- sample(c(-1, 1), n, replace = TRUE)
  B <- V + 0.15 * S
  # The observation probability depends on B beyond V (not the easy MAR case).
  M <- rbinom(n, 1, ifelse(S == 1, 0.8, 0.4))
  phi <- (2 / 3) * truth(V + 0.15) + (1 / 3) * truth(V - 0.15)
  V[M == 0] <- phi[M == 0] <- NA_real_
  list(B = B, V = V, M = M, phi = phi)
}

results <- weights <- curves <- list()
for (scenario_name in names(scenarios)) {
  scenario <- scenarios[[scenario_name]]
  for (problem in c("bridge", "adjoint")) {
    truth <- scenario[[if (problem == "bridge") "beta" else "lambda"]]
    grid <- seq(if (problem == "bridge") -0.95 else -0.90,
                if (problem == "bridge") 0.95 else 0.90, length.out = 301)
    for (n in sizes) for (seed in seeds) {
      dat <- if (problem == "bridge") generate_bridge(n, truth, seed) else generate_adjoint(n, truth, seed)
      timer <- proc.time()[["elapsed"]]
      control <- list(scale = FALSE, bandwidth = 0.30, n_centers = 81L, block_size = 2048L)
      fit <- if (problem == "bridge") {
        fit_bridge_ensemble(dat$B, dat$V, dat$M, scenario$library,
          n_folds = n_folds, seed = 8000L + seed, kernel_control = control)
      } else {
        fit_adjoint_ensemble(dat$B, dat$V, dat$M, dat$phi, scenario$library,
          n_folds = n_folds, seed = 8000L + seed, kernel_control = control)
      }
      elapsed <- proc.time()[["elapsed"]] - timer
      estimates <- vapply(fit$candidates, function(f) predict(f, grid), numeric(length(grid)))
      candidate_rmse <- sqrt(colMeans((estimates - truth(grid))^2))
      ensemble_estimate <- predict(fit, grid)
      correct_index <- match(scenario$correct, names(fit$weights))
      misspecified_rmse <- min(candidate_rmse[-correct_index])
      row <- data.frame(scenario = scenario_name, problem = problem, n = n, seed = seed,
        correct_candidate = scenario$correct,
        selected_candidate = names(fit$weights)[which.max(fit$weights)],
        correct_weight = unname(fit$weights[correct_index]),
        ensemble_rmse = sqrt(mean((ensemble_estimate - truth(grid))^2)),
        correct_candidate_rmse = unname(candidate_rmse[correct_index]),
        best_misspecified_rmse = misspecified_rmse,
        cv_loss = fit$cv_loss, kkt_gap = fit$solver$kkt_gap, elapsed_seconds = elapsed)
      results[[length(results) + 1L]] <- row
      weights[[length(weights) + 1L]] <- data.frame(scenario = scenario_name,
        problem = problem, n = n, seed = seed, candidate = names(fit$weights),
        correctly_specified = names(fit$weights) == scenario$correct,
        weight = unname(fit$weights), cv_loss = unname(fit$candidate_cv_loss),
        rmse = unname(candidate_rmse))
      if (n == max(sizes) && seed == seeds[1L]) {
        curves[[length(curves) + 1L]] <- data.frame(scenario = scenario_name,
          problem = problem, grid = grid, truth = truth(grid), ensemble = ensemble_estimate,
          correct_estimate = estimates[, correct_index],
          misspecified_1 = estimates[, -correct_index, drop = FALSE][, 1L],
          misspecified_2 = estimates[, -correct_index, drop = FALSE][, 2L])
      }
      print(row, row.names = FALSE)
      # Checkpoint every run so failures retain numerical evidence.
      write.csv(do.call(rbind, results), file.path(out_dir, "validation_results.csv"), row.names = FALSE)
      write.csv(do.call(rbind, weights), file.path(out_dir, "candidate_weights.csv"), row.names = FALSE)
    }
  }
}
res <- do.call(rbind, results)
summary <- do.call(rbind, lapply(split(res, list(res$scenario, res$problem, res$n), drop = TRUE), function(x) {
  data.frame(scenario = x$scenario[1L], problem = x$problem[1L], n = x$n[1L],
    seeds = nrow(x), correct_candidate = x$correct_candidate[1L],
    mean_correct_weight = mean(x$correct_weight), min_correct_weight = min(x$correct_weight),
    selection_rate = mean(x$selected_candidate == x$correct_candidate),
    mean_ensemble_rmse = mean(x$ensemble_rmse), max_ensemble_rmse = max(x$ensemble_rmse),
    mean_correct_rmse = mean(x$correct_candidate_rmse),
    min_misspecified_rmse = min(x$best_misspecified_rmse))
}))
summary <- summary[order(summary$scenario, summary$problem, summary$n), ]
write.csv(summary, file.path(out_dir, "validation_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, curves), file.path(out_dir, "seed1_curves.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
saveRDS(list(sizes = sizes, seeds = seeds, n_folds = n_folds,
  libraries = lapply(scenarios, function(x) x$library),
  correct_candidates = vapply(scenarios, function(x) x$correct, character(1)),
  truth_parameters = list(cubic_bridge = c(3, 0.6, 0.6, 1.1),
    cubic_adjoint = c(1.2, 0.4, 0.6, 1.1), rbf_centers = rbf_centers,
    rbf_bandwidth = rbf_bandwidth, rbf_coef = rbf_coef),
  bridge_strength = 0.97, adjoint_delta = 0.15,
  adjoint_observation_probabilities = c(minus = 0.4, plus = 0.8),
  scoring_kernel = control, package_version = as.character(packageVersion("cmbridge"))),
  file.path(out_dir, "configuration.rds"))

pdf(file.path(out_dir, "ensemble_validation.pdf"), width = 10, height = 5)
for (scenario_name in names(scenarios)) for (problem in c("bridge", "adjoint")) {
  ss <- summary[summary$scenario == scenario_name & summary$problem == problem, ]
  par(mfrow = c(1, 2))
  plot(ss$n, ss$mean_correct_weight, type = "b", log = "x", ylim = c(0, 1),
    xlab = "Sample size", ylab = "Mean weight on correct candidate",
    main = paste(scenario_name, problem))
  abline(h = 0.95, lty = 2, col = "grey50")
  curve <- curves[[which(vapply(curves, function(x)
    x$scenario[1L] == scenario_name && x$problem[1L] == problem, logical(1)))]]
  plot(curve$grid, curve$truth, type = "l", lwd = 2,
    xlab = "Function argument", ylab = "Function value", main = "Largest sample, first seed")
  lines(curve$grid, curve$ensemble, col = "blue", lwd = 2, lty = 2)
  legend("topright", c("Truth", "Ensemble"), col = c("black", "blue"), lty = c(1, 2), bty = "n")
}
dev.off()
script_file <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())])
source(file.path(dirname(script_file), "plot_ensemble_truth.R"))
plot_ensemble_truth(out_dir)
print(summary, row.names = FALSE)

# Prespecified acceptance criteria are applied to every largest-sample seed.
# The smaller sizes describe finite-sample behavior; they need not select a vertex.
large <- res[res$n == max(sizes), ]
passed <- all(large$correct_weight >= 0.90) &&
  all(large$selected_candidate == large$correct_candidate) &&
  all(large$ensemble_rmse < 0.06) &&
  all(large$ensemble_rmse < large$best_misspecified_rmse / 2) &&
  all(large$kkt_gap <= 1e-8)
writeLines(c(paste0("passed=", passed), "minimum_correct_weight=0.90",
  "maximum_ensemble_rmse=0.06", "ensemble_rmse_less_than_half_best_misspecified=TRUE"),
  file.path(out_dir, "acceptance.txt"))
if (!passed) stop("Ensemble validation failed its prespecified acceptance criteria.", call. = FALSE)
cat("Ensemble validation passed. Outputs:", out_dir, "\n")
