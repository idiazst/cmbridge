\name{fit_bridge_ensemble}
\alias{fit_bridge_ensemble}
\title{Cross-validated convex ensemble for a bridge equation}
\usage{
fit_bridge_ensemble(B, V, M, library = NULL, n_folds = 5L,
  fold_id = NULL, seed = 1L, kernel = c("rbf", "linear", "cell"),
  kernel_control = list(), solver_control = list(), statistic = c("u", "v"),
  selection = c("rkhs", "robins", "syrgkanis"), selection_control = list())
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
\item{statistic}{\code{"u"} (default) excludes self-products; \code{"v"} retains them. Applies to ensemble validation, candidate penalty CV, and training/new-data moment losses.}
\item{selection}{Default RKHS stacking, or experimental \code{"robins"} or \code{"syrgkanis"} selection. The latter tune all three native methods, including PMMR ridge and Landweber stopping iterations.}
\item{selection_control}{Controls for adaptive selection: \code{projection_basis} (cell_linear), \code{projection_degree} (1), \code{projection_df} (6), positive \code{projection_ridges} (15 log-spaced values from 1e-6 to 10), \code{critic_ridge} (1e-6), \code{bootstrap_draws} (199), \code{noise_probability} (0.95), \code{iteration_grid} (21 log-spaced integer values from 2 to 2000), \code{max_iterations} (100000), \code{max_extensions} (16), \code{extension_factor} (10), \code{extension_points} (4), and \code{limit_tolerance} (1e-8).}
}
\description{Combines candidates to solve \eqn{E[M\beta(V)-1\mid B]=0}, using the common RKHS maximum-moment validation criterion of Section 5.3.2 of the methodology paper.}
\details{The experimental adapters retain the native learner objectives and links.
They are selection principles, not complete reimplementations of the cited
papers' regularized estimators; their theoretical guarantees are not asserted.
Robins scoring uses separate person samples for candidate fitting, conditional
projection training and validation. Its score is twice the product of the
estimated conditional residual and the observed residual, minus the squared
estimated conditional residual, omitting a common constant. Projection ridge
is selected by ordinary leave-one-out prediction loss, with an extendable
positive grid. The known bridge right-hand side is exactly one. Adjoint
conditional loading is estimated in the projection sample. Raw indefinite
Robins Grams are minimized by face enumeration for at most ten candidates;
they are not projected to PSD.

Syrgkanis penalty selection uses the actual full-training fit and its
covariance-stabilized empirical critic. Its threshold is calibrated for that
sample size, rather than transported from smaller internal fits. A common simultaneous Gaussian
multiplier-bootstrap envelope over the entire evaluated path calibrates the
noise threshold. The largest admissible regularization is selected. An
unbracketed threshold extends the grid. A converged most-regularized limit,
or Landweber's natural one-iteration endpoint, can be accepted and is recorded
explicitly; this is not acceptance of a boundary CV minimum. Ensemble weights
minimize the empirical second moment of the fitted function subject to the
held-out discrepancy threshold. No feasible mixture raises an error.

Adaptive selections use their projection basis instead of kernel/statistic
for selection. Their training and new-data \code{moment_loss} remains a
separate RBF U-statistic diagnostic. Supplied kernel/statistic parameters only
control the default RKHS selection pathway. All projection, basis, threshold
and penalty fitting remains inside the supplied training sample. Return
values record disjoint selection samples, projection tuning, grids, threshold
events and ensemble constraint status. The adjoint loading callback is also
used for adaptive training splits.

Every kernel supports U and V scoring: cell, Gaussian, and linear.
Penalty CV uses the same statistic, kernel and training-only kernel specification as
ensemble selection. Failed penalty trials are recorded and never selected.
Candidates which fail tuning, fitting, or validation are excluded across all
folds and given weight zero; the returned object retains their failure details.
An ensemble with no valid candidates raises an error.

For each inner fold, fit candidates on the other folds and form held-out residuals. Gaussian and linear V scores use \eqn{G_f=R_f^T K_f R_f/n_f^2}; U scores replace the kernel diagonal by zero and divide by \eqn{n_f(n_f-1)}. When \code{M=0}, the residual is -1 and V is not evaluated. The common Gram matrix is the average of fold-specific matrices weighted by their scored sample sizes. Weights minimize \eqn{\alpha^T G\alpha} over nonnegative weights summing to one. Candidates are refit on the entire supplied sample before prediction.

Nystr\"om scoring uses a common low-rank approximation of the Gaussian kernel, independently of the kernels or bases used to fit the candidate learners. Exact scoring evaluates the full Gaussian kernel with the selected statistic. The linear kernel is a finite-feature conditional-moment diagnostic, not a generally identifying kernel.

To preserve outer cross-fitting, call this function with only the outer training observations and predict on the untouched outer validation sample. This function performs inner stacking; it does not construct the complete longitudinal causal estimator.

Cell U scoring uses distinct pairs within each joint conditioning cell. For a
cell containing \eqn{m} observations in a validation sample of size \eqn{n},
its contribution is the sum of ordered residual cross-products with unequal
indices, divided by \eqn{n(m-1)}. Singleton cells contribute zero, including samples where every
cell is a singleton. Cell V scoring includes all ordered pairs and divides by
\eqn{nm}; singleton cells therefore contribute their residual outer product
divided by \eqn{n}. All conditioning columns are retained.
Raw U-statistic matrices can be indefinite. The averaged Gram is projected
onto the positive semidefinite cone by replacing negative eigenvalues by zero
before simplex optimization. V matrices are already positive semidefinite and
are used directly without projection. The raw matrix, raw candidate losses,
projection adjustment and whether projection was applied are retained. Penalty
selection uses the raw scalar selected statistic, without truncating negative
U losses to zero. V self-products introduce a candidate-dependent finite-sample
bias; nonnegative scores do not imply more accurate weight selection.

A sieve, Landweber or saturated L1 candidate can supply \code{penalty_scales}. For
\code{sieve_md}, these are absolute ridge \code{lambda} values. For
\code{landweber}, these are absolute positive \code{weight_ridge} values
regularizing the conditioning-moment weight matrix, not the target coefficients.
For
\code{saturated_l1}, each scale is divided by the square root of the number
of training people to obtain \code{penalty}. Penalty selection
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

Each refitted candidate retains its selected \code{penalty_cv},
\code{penalty_boundary}, \code{penalty_status}, \code{penalty_parameter},
\code{penalty_ids}, and \code{penalty_fold_id}. The controls
\code{penalty_ids} and \code{penalty_folds(ids)} let an outer estimator pass
the same person splits used by its other learners. Penalties are chosen
separately inside every ensemble training fold and again for the final refit.

}
\value{A \code{cmbridge_ensemble} inheriting from \code{cmbridge_fit}. Fields include \code{raw_gram}, \code{psd_projection}, \code{candidate_raw_cv_loss}, \code{raw_cv_loss}, \code{fold_penalty_cv}, \code{fold_penalty_boundary}, \code{weights}, refitted \code{candidates}, \code{library}, \code{fold_id}, \code{gram}, \code{fold_gram}, \code{fold_n_scored}, \code{cv_predictions}, \code{cv_residuals}, \code{candidate_cv_loss}, \code{cv_loss}, and solver diagnostics. Prediction and \code{moment_loss} are supported. Cross-validated predictions on unobserved V are NA. The stored training moment loss uses a kernel learned from the full supplied training sample.}
\seealso{\code{\link{fit_bridge}}, \code{\link{fit_adjoint_ensemble}}, \code{\link{moment_loss}}}
