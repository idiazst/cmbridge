\name{fit_cm}
\alias{fit_cm}
\title{Fit a conditional-moment function}
\usage{fit_cm(response, diagonal = NULL, target, instrument,
  method = c("sieve_md", "landweber", "pmmr", "saturated_l1"), control = list())}
\arguments{
\item{response}{Numeric response vector.}
\item{diagonal}{Numeric multiplicative/diagonal vector. Defaults to one.}
\item{target}{Variables indexing the unknown function.}
\item{instrument}{Variables defining the conditional moment.}
\item{method}{Estimation method.}
\item{control}{Method-specific control list. The default identity link is unbounded. Finite \code{lower} or \code{upper} values raise an error; use the link controls below for a model parameterization inside the optimization. Post-fit clipping is not provided.}
}
\description{Fits h in E[Y-D h(X) given Z]=0 using a common interface.}
\details{For \code{sieve_md} and \code{landweber}, \code{target_basis} and
\code{instrument_basis} may be \code{"bs"}, \code{"poly"}, \code{"quadratic"}, \code{"cell"}, or \code{"cell_linear"}.
The \code{"quadratic"} basis includes an intercept, every main effect, every
square, and every pairwise interaction. These terms are fixed independently
of training-sample combinations, with no category fallback.
The \code{"cell"} basis includes an intercept and indicators for joint input
categories observed in training, representing all their interactions; unseen
categories receive the intercept. When both sieve bases are \code{"cell"},
grouped sufficient statistics evaluate the same normal equations without
constructing dense observation-by-category matrices.
With both sieve bases set to \code{"cell"}, \code{link = "inverse_logit"}
optimizes unrestricted real cell coefficients and predicts
\eqn{1/\operatorname{expit}(\eta) = 1+\exp(-\eta)}.
This gives values strictly above one mathematically, with one as a limit.
\code{link = "log"} uses \eqn{\exp(\eta)} for positivity.
The default \code{link = "identity"} is unrestricted on the function scale.
The nonlinear solver minimizes the original squared conditional-moment
objective plus its ridge penalty on function-scale cell deviations from an
unpenalized common value. It uses unconstrained BFGS, with unconstrained
\code{nlminb} refinement when needed. It imposes no coefficient bounds and
never clips a fitted function. Unseen cells receive the linked common value.
Controls are \code{solver_tolerance} (1e-10) and \code{solver_max_iter} (10000).
\code{solver} records optimizer traces, gradients, and convergence failures.
The optimized unrestricted coefficients are stored in \code{link_coefficients}
and \code{solver$parameters}. The \code{coefficients} field retains the
function-scale common value and cell deviations for existing diagnostics.
Positive sieve links also support the other basis choices, using ridge on scaled
link coefficients and an unpenalized intercept. The \code{"cell_linear"} basis
has all fixed main effects followed by joint-category deviations, which receive
\code{cell_penalty_multiplier} times the ridge (default 100). Unseen categories
retain the main-effect prediction. Pure cell fits may request
\code{penalty_on = "link"}; their default remains the function-value penalty.
Large function-value penalties retain the common direction by a separate solve.
For \code{landweber}, the same links may be used with any supported basis;
unrestricted nonlinear gradient updates with a fixed initial step and backtracking minimize the
conditional-moment objective, regularized by iteration stopping.
\code{precondition = TRUE} makes an invertible change of coordinates using the
training design covariance; this preserves every term in the specified class.
For linked Landweber splines, \code{target_extrapolation} defaults to
\code{"constant"}: the basis continues at its boundary outside the training
range. This is part of the model, not clipping of fitted function values.
Identity fits retain \code{"polynomial"} continuation.
For \code{pmmr}, the links transform the target Nystrom predictor inside
the moment objective. The RKHS penalty applies to that linear predictor,
rather than to the transformed function. Its unrestricted BFGS solver
records optimizer traces and checks the coefficient gradient.
These links are opt-in through \code{control$link}; the generic defaults
remain identity for compatibility. The lmtp bridge integration selects
inverse-expit for all three bridge candidates and identity for adjoints.
The \code{saturated_l1} method expands joint target and instrument cells
and minimizes weighted squared empirical conditional moments with an L1
penalty on target-cell deviations from an unpenalized constant function.
The constant is multiplied by D. Controls are \code{penalty} (0.001),
\code{target_columns} and \code{instrument_columns} (all columns),
\code{lower}/\code{upper} (infinite), \code{tolerance} (1e-8), and
\code{max_iter} (1000000). Unseen target cells use the fitted constant.
This categorical learner assumes fixed finite support; numeric values define
distinct joint cells. Solver failure is an error.}
\value{An object of class \code{cmbridge_fit}.}
