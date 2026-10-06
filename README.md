# cmbridge

`cmbridge` standardizes learners for conditional moment equations of the form

`E[Y - D h(X) | Z] = 0`.

The initial learner library contains:

1. `sieve_md`: sieve minimum distance with polynomial, additive B-spline, or saturated joint-category bases;
2. `landweber`: Landweber iterative regularization of the same empirical conditional-moment operator;
3. `pmmr`: scalable Gaussian-kernel proxy/maximum-moment restriction using low-rank target and critic features.

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
