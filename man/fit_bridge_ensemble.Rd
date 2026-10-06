\name{fit_bridge_ensemble}
\alias{fit_bridge_ensemble}
\title{Cross-validated convex ensemble for a bridge equation}
\usage{
fit_bridge_ensemble(B, V, M, library = NULL, n_folds = 5L,
  fold_id = NULL, seed = 1L, kernel = c("rbf", "linear", "cell"),
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
\item{kernel}{Common validation kernel: Gaussian RBF, \code{linear} (an intercept plus standardized linear features), or \code{cell} (joint categorical cells).}
\item{kernel_control}{List with \code{scale} (default TRUE), \code{bandwidth} (default training-sample median distance), \code{approximation} (\code{"nystrom"}, default, or \code{"exact"}), \code{n_centers} (100), \code{block_size} (512), optional \code{seed}, \code{nystrom_ridge} (1e-10), and \code{nystrom_tol} (1e-10). Gaussian centers, scaling, and bandwidth are learned exclusively in each inner training sample. Increase the rank as needed; a fixed low-rank kernel need not detect every conditional-moment violation. Exact Gaussian scoring uses blocks to avoid an n-by-n stored kernel, but retains quadratic computation.}
\item{solver_control}{List with \code{tol} (1e-8) and \code{max_iter} (50000) for simplex quadratic optimization. The tolerance applies after scaling the Gram matrix by its largest absolute entry.}
}
\description{Combines candidates to solve \eqn{E[M\beta(V)-1\mid B]=0}, using the common RKHS maximum-moment validation criterion of Section 5.3.2 of the methodology paper.}
\details{
For each inner fold, fit candidates on the other folds, form held-out residuals, and compute \eqn{G_f=R_f^T K_f R_f/n_f^2}. When \code{M=0}, the residual is -1 and V is not evaluated. The common Gram matrix is the average of fold-specific matrices weighted by their scored sample sizes. Weights minimize \eqn{\alpha^T G\alpha} over nonnegative weights summing to one. Candidates are refit on the entire supplied sample before prediction.

Nystr\"om scoring uses a common low-rank approximation of the Gaussian kernel, independently of the kernels or bases used to fit the candidate learners. Exact scoring evaluates the full Gaussian V-statistic. The linear kernel is a finite-feature conditional-moment diagnostic, not a generally identifying kernel.

To preserve outer cross-fitting, call this function with only the outer training observations and predict on the untouched outer validation sample. This function performs inner stacking; it does not construct the complete longitudinal causal estimator.

Cell scoring uses distinct pairs within each joint conditioning cell. For a
cell containing \eqn{m} observations in a validation sample of size \eqn{n},
its contribution is the sum of ordered residual cross-products with unequal
indices, divided by \eqn{n(m-1)}. Singleton cells contribute zero; if every
cell is a singleton, scoring fails. All conditioning columns are retained.
Raw U-statistic matrices can be indefinite. The averaged Gram is projected
onto the positive semidefinite cone by replacing negative eigenvalues by zero
before simplex optimization. The raw matrix, raw candidate losses, and the
projection adjustment are retained. Penalty selection uses the raw scalar
U-statistic, without truncating negative losses to zero.

A saturated L1 candidate can supply \code{penalty_scales}. Penalty selection
uses person-based training-only splits, and records boundary minima as failed
attempts. It extends the grid by one decade with four intermediate logarithmic
points on the indicated side until a minimum is interior. The controls
\code{penalty_max_extensions} (16), \code{penalty_extension_factor} (10), and
\code{penalty_extension_points} (4) configure this process. An unresolved
boundary raises \code{cmbridge_penalty_boundary_error} with the complete
\code{penalty_cv} and \code{boundary_history}. An exact tie is accepted only
when an interior grid point also attains the minimum. Include
\code{method="saturated_l1"} explicitly for discrete bridge problems; the
original three-method default is retained.

}
\value{A \code{cmbridge_ensemble} inheriting from \code{cmbridge_fit}. Fields include \code{raw_gram}, \code{psd_projection}, \code{candidate_raw_cv_loss}, \code{raw_cv_loss}, \code{fold_penalty_cv}, \code{fold_penalty_boundary}, \code{weights}, refitted \code{candidates}, \code{library}, \code{fold_id}, \code{gram}, \code{fold_gram}, \code{fold_n_scored}, \code{cv_predictions}, \code{cv_residuals}, \code{candidate_cv_loss}, \code{cv_loss}, and solver diagnostics. Prediction and \code{moment_loss} are supported. Cross-validated predictions on unobserved V are NA. The stored training moment loss uses a kernel learned from the full supplied training sample.}
\seealso{\code{\link{fit_bridge}}, \code{\link{fit_adjoint_ensemble}}, \code{\link{moment_loss}}}
