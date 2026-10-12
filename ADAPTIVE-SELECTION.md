# Experimental conditional-moment selection

Branch `robins-syrgkanis-selection` adds `selection="robins"` and
`selection="syrgkanis"` to `fit_cm`, `fit_bridge`, `fit_adjoint`, and both
ensemble interfaces. Defaults remain fixed controls for individual learners
and the existing RKHS U/V selector for ensembles.

Both adapters tune sieve ridge, PMMR ridge, and Landweber's actual stopping
iterations. Positive ridge grids are used throughout. Inverse-expit bridge
parameterizations and unrestricted adjoints are unchanged. No predictions
are clipped, and no outcome-regression or TMLE update is changed.

## Robins

Adapt Eq. 6 of [Ghassami, Robins and Rotnitzky (2025)](https://arxiv.org/abs/2505.20787).
Each split has three person-disjoint samples: estimator fitting, conditional
projection training, and scoring. Supplied two-fold inner plans are respected:
one side fits the candidate, the other is subdivided into projection/scoring
samples, with both orientations used. Three or more folds rotate these roles.
Projection models use a saturated joint-cell plus main-effect ridge design,
with training-only, positive, extendable leave-one-out ridge selection.

For a candidate output U = D*h(X), conditional projection t, loading Y and
conditional loading projection r, the score is
`mean(2*(t-r)*(U-Y) - (t-r)^2)`. This equals the paper's score minus the common
constant `mean((r-Y)^2)`. For the bridge Y=r=1 is known. Adjoint loading
projections are estimated. Mixtures combine candidate projections linearly.
Raw score matrices are optimized globally by simplex face enumeration for
libraries of at most ten candidates, without PSD projection.

## Syrgkanis

Adapt the discrepancy principle of [Tan and Syrgkanis (2026)](https://arxiv.org/abs/2603.01337).
Native fits retain their native objectives and penalties. Training equation
error is evaluated on the actual full-training fit, with a critic basis
and second-moment matrix learned on those same training rows, plus a positive
ridge. It is a training discrepancy rule, not penalty cross-validation.
The threshold is a simultaneous multiplier-bootstrap envelope of centered
moment fluctuations across every evaluated penalty, using common bootstrap
draws for all penalties within a split. It is estimated from data, never DGP
truth. The strongest admissible regularization is chosen; a missing crossing
extends the path. A converged most-regularized limit is reported explicitly.
The selected penalty is calibrated for that full-training sample size, rather
than transported from smaller internal fits. This endpoint is different from a minimum at a CV grid boundary: Robins and
projection CV still require interior minima.

The paper does not specify our ensemble problem. Its adaptation minimizes
the empirical second moment of the candidate mixture subject to a held-out
discrepancy threshold, using a regularization-multiplier path. If no mixture
satisfies that threshold, selection fails with its minimum loss and threshold
preserved. Thresholds are not increased to make failed fits pass.

## Limits

These are adapters around the existing native learners, not implementations
of the papers' complete regularized estimators. In particular, the papers'
function-scale penalties and source-condition assumptions do not automatically
cover our link-scale penalties. The bootstrap calibration is an explicit
practical adaptation, not the papers' class-complexity bound. Neither score
equals full function MSE, nor resolves nonuniqueness.

All final candidate fits use the entire supplied training sample. The
independence needed to learn projections reduces the sample available for
intermediate candidate fits relative to the old two-way RKHS stacking.
Comparisons therefore change the split roles and tune previously fixed
controls, as well as change the selection criterion.

Training/new-data `moment_loss()` remains an RBF U diagnostic, separate from
the adaptive selection score. Inspect `selection_audit`, `fold_scores`,
`penalty_cv`, `penalty_boundary`, and `solver` for the actual criterion and
calibration details. No complete causal-estimator coverage claim follows from
the limited bridge/adjoint diagnostic.
