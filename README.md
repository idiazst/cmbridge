# cmbridge

`cmbridge` standardizes three learners for conditional moment equations of the form

`E[Y - D h(X) | Z] = 0`.

The initial learner library contains:

1. `sieve_md`: sieve minimum distance with polynomial or additive B-spline bases;
2. `landweber`: Landweber iterative regularization of the same empirical conditional-moment operator;
3. `pmmr`: scalable Gaussian-kernel proxy/maximum-moment restriction using low-rank target and critic features.

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

See `REVIEW.md` for the implementation review and method-selection rationale, `REFERENCES.bib` for citations, and `inst/simulations/validate_large_sample.R` for the large-sample recovery study.
