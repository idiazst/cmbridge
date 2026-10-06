# cmbridge

`cmbridge` standardizes learners for conditional moment equations of the form

`E[Y - D h(X) | Z] = 0`.

The initial learner library contains:

1. `sieve_md`: sieve minimum distance with polynomial, additive B-spline, or saturated joint-category bases;
2. `landweber`: Landweber iterative regularization of the same empirical conditional-moment operator;
3. `pmmr`: scalable Gaussian-kernel proxy/maximum-moment restriction using low-rank target and critic features.

Landweber also accepts `target_basis = "quadratic"`: an intercept, all main
effects, squares, and every pairwise interaction of the supplied variables.
Its terms are fixed independently of the observed combinations. With
`link = "inverse_logit"`, these terms parameterize the linear predictor of
the bridge; the optimization coefficients remain unrestricted.

For finite discrete support, `saturated_l1` adds all joint-cell interactions and
fits the empirical conditional moments with an L1 penalty. The original three
learners remain the default ensemble library.

The package exposes a common generic interface and two manuscript-oriented wrappers:

```r
# Generic conditional moment
fit <- fit_cm(response = Y, diagonal = D,
              target = X, instrument = Z,
              method = "sieve_md")

# Inverse-probability bridge: E[M beta(V)-1 | B] = 0
fit_beta <- fit_bridge(B = B, V = V, M = M, method = "pmmr")

# Adjoint: E[M{lambda(B)-phi(V)} | V] = 0
fit_lambda <- fit_adjoint(B = B, V = V, M = M, phi = phi,
                          method = "landweber")

predict(fit_beta, newdata)
moment_loss(fit_beta)
```

`V` is allowed to be missing when `M=0` in `fit_bridge()`. The adjoint wrapper uses complete cases, which is equivalent because its conditional moment is multiplied by `M`.

## Installation

```r
install.packages("remotes")
remotes::install_github("idiazst/cmbridge")
library(cmbridge)
```

## Cross-validated ensembles

Use `control = list(target_basis = "cell", instrument_basis = "cell")`
with `sieve_md` or `landweber` to give each observed joint combination of
discrete inputs its own coefficient, including all interactions. The basis
uses an unpenalized common constant plus cell deviations; previously unseen
target combinations receive that constant. Raw predictions are unbounded. Finite `lower` or `upper` controls are rejected;
post-fit clipping can invalidate a conditional-moment solution.

The ensemble functions implement the convex conditional-moment stacking in
Section 5.3.2 of the methodology paper. Fit candidate specifications within
inner training folds, score their held-out residuals with one common RKHS
kernel, solve the simplex quadratic program, and refit all candidates on the
entire supplied training sample:

```r
library <- list(
  cubic = list(method = "sieve_md", control = list(
    target_basis = "poly", target_degree = 3L)),
  iterative = list(method = "landweber"),
  kernel = list(method = "pmmr")
)
ensemble <- fit_bridge_ensemble(B, V, M, library = library, n_folds = 5L)
ensemble$weights
predict(ensemble, newdata = V_new)

# phi_known is a known or independently estimated loading vector.
adjoint <- fit_adjoint_ensemble(B, V, M, phi_known, library = library)
predict(adjoint, newdata = B_new)
```

For a loading estimated from the supplied sample, pass a callback
`phi(train, validation)` returning `list(train = ..., validation = ...)`.
Estimate the loading using only the supplied training indices and evaluate
that same loading at both sets of indices. The callback is called once per
inner fold and once for the final refit (with empty validation indices).
All candidates in a fold receive the same loading. Adjoint scoring uses only
validation complete cases; bridge scoring retains all rows, assigning residual
`-1` when `M=0` without evaluating missing `V`.

By default the common scorer uses a training-derived Nyström approximation
of a Gaussian kernel, with 100 centers, to support large samples. It is
independent of candidate fitting features. Increase `n_centers` as needed;
a fixed low-rank kernel need not detect all conditional-moment violations.
For the exact Gaussian V-statistic, use
`kernel_control = list(approximation = "exact")`; its memory use is blocked
but its computational cost is quadratic in the scored sample size.
Fold-specific Gram matrices are averaged with weights proportional to their
scored sample sizes. `cv_loss`, `candidate_cv_loss`, and `gram` expose the
criterion used to select the weights; `moment_loss()` supports the common
scorer on new data.

When using an ensemble in an outer cross-fitted causal estimator, supply
only that outer training sample to the ensemble function, then predict on
the untouched outer validation sample. These wrappers perform inner
stacking and do not implement the full longitudinal causal estimator.

The independent ensemble selection validation is in
[`idiazst/cmbridge-tests`](https://github.com/idiazst/cmbridge-tests).

## Saturated discrete learners

```r
cells <- list(
  small = list(method = "saturated_l1", control = list(penalty = 0.001 / sqrt(length(M)))),
  medium = list(method = "saturated_l1", control = list(penalty = 0.01 / sqrt(length(M)))),
  large = list(method = "saturated_l1", control = list(penalty = 0.1 / sqrt(length(M))))
)
bridge <- fit_bridge_ensemble(B, V, M, library = cells,
                              fold_id = shared_inner_labels, kernel = "cell")
```

The cell scorer uses ordered pairs of distinct observations within each joint
conditioning cell. It excludes self-products and divides each cell's sum by
`n * (cell_count - 1)`. Singleton cells contribute zero and their count is
recorded; a validation sample containing only singleton cells raises an error.
Every conditioning variable, including treatment, is retained. Ensemble weights
use the PSD projection of the averaged U-statistic Gram. `raw_gram`,
`candidate_raw_cv_loss`, `raw_cv_loss`, and `psd_projection` retain the original
criterion and projection diagnostics. Penalty CV uses raw U-statistic losses,
which can be negative.

Positive penalty grids are extended when their minimum is at a boundary. Every
boundary choice is recorded as a failed attempt. Only an interior minimum is
accepted; an unresolved boundary raises an error carrying the losses and attempt
history. The default limit is 16 one-decade extensions, with four logarithmic
points per extension. No zero penalty is introduced. `select_penalty_grid()` is
shared with lmtp, and `cell_moment_gram()` supplies the common scorer.

For a discrete `sieve_md` ensemble candidate, supply
`control$penalty_scales = 10^seq(-10, 0, by = .25)` to select the ridge
`lambda` on training-only CV splits. These are absolute lambda values.
`penalty_ids` and the `penalty_folds(ids)` callback preserve shared person
splits, including when tuning inside an ensemble training fold. The selected
grid, extension attempts, person IDs, and CV labels are saved in each fit.
Adjoint tuning reconstructs its loading within each penalty training sample.
The same controls continue to support `saturated_l1`, whose scales are divided
by the square root of the number of training people.

For a `landweber` candidate, `penalty_scales` selects `weight_ridge`, the
positive ridge in the conditioning-moment weight matrix. For example,
`penalty_scales = 10^seq(-6, 0, by = .5)` supplies thirteen initial values.
This regularizes the relative weights of conditioning directions; it does
not penalize the target coefficients or remove conditioning variables.
The same training-only person splits, raw U-statistic validation scores,
boundary extension and failure rules apply. `penalty_parameter` records
`"weight_ridge"`, and no tuning occurs unless the grid is explicitly supplied.

Learners use an
unpenalized constant function and penalize joint-cell deviations; the constant
is multiplied by the diagonal in the conditional-moment equation.

Controls include `target_columns`, `instrument_columns`, `tolerance`, and
`max_iter`. Predictions are unbounded; unseen target cells use the fitted
constant. Finite bounds require a constrained optimizer and are rejected here. A warm-start
penalty path is used, and incomplete solver fits raise an error.

The observed-history longitudinal integration and exact-truth simulation are
in the sibling lmtp development package and cmbridge-tests repository. Broad
function classes contain the true discrete functions; small effective cell
counts and weak inverse operators can still affect finite-sample performance.

See `REVIEW.md` for the implementation review and method-selection rationale, `REFERENCES.bib` for citations, and `inst/simulations/validate_large_sample.R` for the large-sample recovery study.

The joint-category `sieve_md` supports `control = list(target_basis = "cell", instrument_basis = "cell", link = "inverse_logit")`. It directly optimizes unrestricted real coefficients and predicts `1 / expit(eta) = 1 + exp(-eta)`, giving values above one. `link = "log"` uses an exponential parameterization, and the default `identity` link remains unrestricted on the function scale. The linked fit keeps the original ridge penalty on function-scale cell deviations, records optimizer traces and gradient checks, and imposes no coefficient bounds or prediction clipping.


Landweber and PMMR also support `link = "inverse_logit"`. Landweber uses
unrestricted nonlinear gradient updates with backtracking and iteration
stopping. PMMR minimizes the transformed empirical moment loss and penalizes
the RKHS linear predictor; it saves optimizer traces and a gradient check.
Generic `fit_cm`/`fit_bridge` defaults retain identity for compatibility with
existing polynomial and kernel specifications. The longitudinal integration
explicitly uses inverse-expit for all three bridge candidates and identity
for all three adjoint candidates.

Linked Landweber splines use constant continuation outside their training
boundaries (`target_extrapolation = "constant"`), and a fixed initial step
that may decrease through backtracking. A regression test reproduces a small
nested numerical-dose split where polynomial spline continuation produced
validation values around 1e189 after the exponential link. The new model
continuation prevents that overflow without clipping fitted bridge values.
Identity adjoint fits retain their previous spline continuation.

CV records include scored and paired row counts and an `informative` flag.
Singleton conditioning cells contribute zero, including a validation subset
with no observed pair. When all scores tie, the documented interior tie rule
is used; an uninformative split does not establish that its penalty is optimal.

The ensemble now removes self-products for Gaussian and linear kernels as well
as cell kernels, and penalty CV uses the same kernel as ensemble selection.
Raw U-statistics can be negative; PSD projection is applied only before the
convex weight optimization.

Linked sieve fits also support fixed `poly`/`quadratic` bases and `cell_linear`
(fixed main effects plus joint-category deviations). Their ridge acts on scaled
link coefficients, leaves the intercept unpenalized, and uses
`cell_penalty_multiplier = 100` for the joint deviations in `cell_linear`.
Every predictor is retained. Pure cell fits can request `penalty_on = "link"`;
the default pure-cell function penalty remains available. An unrestricted
damped Gauss-Newton solver and BFGS refinement record gradient diagnostics.

`landweber` supports `precondition = TRUE`, an invertible change of coordinates
which retains the entire specified basis class. Cell critic redundancy is
removed by an exact ridge-weight equivalence. Common cell directions are solved
separately from large ridge penalties to retain numerical precision.

Failed penalty trials (including nonfinite scores) are retained in CV records
and never selected. A candidate that cannot pass tuning/fitting/validation is
recorded in `candidate_failures` and excluded across all folds and the final
refit, with weight zero. All other candidates retain the same validation rows.
An ensemble with no valid candidates raises a structured error. Boundary minima
still extend the grid or fail; this does not accept a boundary minimum.
