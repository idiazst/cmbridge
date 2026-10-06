\name{fit_adjoint_ensemble}
\alias{fit_adjoint_ensemble}
\title{Cross-validated convex ensemble for an adjoint equation}
\usage{
fit_adjoint_ensemble(B, V, M, phi, library = NULL, n_folds = 5L,
  fold_id = NULL, seed = 1L, kernel = c("rbf", "linear", "cell"),
  kernel_control = list(), solver_control = list())
}
\arguments{
\item{B}{Numeric vector or matrix indexing the adjoint function.}
\item{V}{Numeric vector or matrix conditioning the adjoint moment. May be missing when \code{M=0}.}
\item{M}{Binary response/visit indicator.}
\item{phi}{A numeric loading vector of length n, or a function \code{phi(train, validation)} returning a list with numeric \code{train} and \code{validation} vectors matching those observation-index vectors. A callback must estimate the loading using only \code{train} and use that same estimate for every candidate and both returned vectors. It is called once per fold and once for the final refit, with all training indices and \code{integer(0)} validation indices. Numeric loadings should be known or estimated independently of the supplied sample; use a callback for loadings learned within this sample. Loadings may be NA on rows with M=0.}
\item{library}{Candidate library, as in \code{fit_bridge_ensemble}, with controls passed to \code{fit_adjoint}.}
\item{n_folds}{Number of inner folds, stratified by M.}
\item{fold_id}{Optional integer fold labels per input observation.}
\item{seed}{Nonnegative integer random seed. The calling random state is restored.}
\item{kernel}{Common scoring kernel on V.}
\item{kernel_control}{Scoring-kernel controls; see \code{fit_bridge_ensemble}.}
\item{solver_control}{Simplex optimization controls; see \code{fit_bridge_ensemble}.}
}
\description{Fits a convex ensemble for \eqn{E[M\{\lambda(B)-\phi(V)\}\mid V]=0}.}
\details{The equivalent complete-case conditional moment is scored only on validation rows with M=1, using residuals \eqn{\lambda(B)-\phi(V)}. Every candidate in a fold receives the same training loading; the loading callback allows its fit to exclude that fold. Fold Gram matrices are weighted by their complete-case counts. Candidates are finally refit with one common full-training-sample loading. All fitting and weight selection must take place inside an outer training sample when this ensemble is used in a cross-fitted causal estimator.}
\value{A \code{cmbridge_ensemble}; see \code{fit_bridge_ensemble}. Cross-validated predictions and residuals are NA when M=0. \code{complete_case_index} records the scored observations.}
\seealso{\code{\link{fit_adjoint}}, \code{\link{fit_bridge_ensemble}}}
