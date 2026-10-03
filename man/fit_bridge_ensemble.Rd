\name{fit_bridge_ensemble}
\alias{fit_bridge_ensemble}
\title{Cross-validated convex ensemble for a bridge equation}
\usage{
fit_bridge_ensemble(B, V, M, library = NULL, n_folds = 5L,
  fold_id = NULL, seed = 1L, kernel = c("rbf", "linear"),
  kernel_control = list(), solver_control = list())
}
\arguments{
\item{B}{Numeric vector or matrix indexing the bridge conditional moment.}
\item{V}{Numeric vector or matrix indexing the bridge function. May be missing when \code{M=0}.}
\item{M}{Binary response/visit indicator.}
\item{library}{Named list of candidates. Each is a list with \code{method} and optional \code{control}, as accepted by \code{fit_bridge}. Character method names are also accepted. The default contains all three methods with their default controls. Tuning configurations can be separate candidates.}
\item{n_folds}{Number of inner folds, stratified by \code{M}.}
\item{fold_id}{Optional integer fold labels, one per input observation. Overrides \code{n_folds}. Every fold needs complete cases in both its training and validation samples.}
\item{seed}{Nonnegative integer random seed. The calling session's random state is restored.}
\item{kernel}{Common validation kernel: Gaussian RBF or \code{linear} (an intercept plus standardized linear features).}
\item{kernel_control}{List with \code{scale} (default TRUE), \code{bandwidth} (default training-sample median distance), \code{approximation} (\code{"nystrom"}, default, or \code{"exact"}), \code{n_centers} (100), \code{block_size} (512), optional \code{seed}, \code{nystrom_ridge} (1e-10), and \code{nystrom_tol} (1e-10). Gaussian centers, scaling, and bandwidth are learned exclusively in each inner training sample. Increase the rank as needed; a fixed low-rank kernel need not detect every conditional-moment violation. Exact Gaussian scoring uses blocks to avoid an n-by-n stored kernel, but retains quadratic computation.}
\item{solver_control}{List with \code{tol} (1e-8) and \code{max_iter} (50000) for simplex quadratic optimization. The tolerance applies after scaling the Gram matrix by its largest absolute entry.}
}
\description{Combines candidates to solve \eqn{E[M\beta(V)-1\mid B]=0}, using the common RKHS maximum-moment validation criterion of Section 5.3.2 of the methodology paper.}
\details{
For each inner fold, fit candidates on the other folds, form held-out residuals, and compute \eqn{G_f=R_f^T K_f R_f/n_f^2}. When \code{M=0}, the residual is -1 and V is not evaluated. The common Gram matrix is the average of fold-specific matrices weighted by their scored sample sizes. Weights minimize \eqn{\alpha^T G\alpha} over nonnegative weights summing to one. Candidates are refit on the entire supplied sample before prediction.

Nystr\"om scoring uses a common low-rank approximation of the Gaussian kernel, independently of the kernels or bases used to fit the candidate learners. Exact scoring evaluates the full Gaussian V-statistic. The linear kernel is a finite-feature conditional-moment diagnostic, not a generally identifying kernel.

To preserve outer cross-fitting, call this function with only the outer training observations and predict on the untouched outer validation sample. This function performs inner stacking; it does not construct the complete longitudinal causal estimator.
}
\value{A \code{cmbridge_ensemble} inheriting from \code{cmbridge_fit}. Fields include \code{weights}, refitted \code{candidates}, \code{library}, \code{fold_id}, \code{gram}, \code{fold_gram}, \code{fold_n_scored}, \code{cv_predictions}, \code{cv_residuals}, \code{candidate_cv_loss}, \code{cv_loss}, and solver diagnostics. Prediction and \code{moment_loss} are supported. Cross-validated predictions on unobserved V are NA. The stored training moment loss uses a kernel learned from the full supplied training sample.}
\seealso{\code{\link{fit_bridge}}, \code{\link{fit_adjoint_ensemble}}, \code{\link{moment_loss}}}
