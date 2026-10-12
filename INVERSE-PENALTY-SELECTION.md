# Inverse-problem parameter selection

The `inverse-penalty-selection` branch reverts the Robins/Syrgkanis adapters.
It adds three optional, training-only parameter-choice rules. Defaults and the
RKHS U/V ensemble objectives, PSD projection and simplex optimization are unchanged.

Use `select_cm_penalty(..., rule="gcv")` directly, or put
`penalty_selection="gcv"` (also `"lcurve"` or `"quasiopt"`) and optional
`penalty_selection_control=list(...)` in a learner's `control`. The latter
works in bridge/adjoint wrappers and each member of an ensemble library.
All three native learners are supported. Sieve and PMMR tune positive ridge;
Landweber tunes its actual positive stopping iteration, not critic ridge.

## Rules

* Generalized cross-validation (GCV) divides squared native equation residual
  by a degrees-of-freedom correction. The equation is expressed in independent
  weighted critic coordinates; redundant columns are removed orthogonally,
  preserving the fitted loss. The trace is exact for linear ridge and fixed
  linear Landweber iterations, including the initialized intercept. For linked
  ridge it uses the local Hessian of the actual penalized objective, including
  residual and penalty curvature. For nonlinear Landweber, central differences
  of the actual iteration include initialization, preconditioning and
  backtracking. Low-dimensional critic traces use all coordinate directions;
  larger traces use common random probes, with the estimated trace standard
  error retained. Invalid traces yield inadmissible GCV trials.
* L-curve selects maximum positive curvature of log equation residual versus
  log native regularization seminorm. The native link-coefficient penalty is
  measured for linked sieve/PMMR; function-centered values are measured for the
  cell sieve. For Landweber the norm is the coefficient deviation from its
  initialization, in the coordinates used by the iteration. Local quadratic
  derivatives on the logarithmic path approximate curvature.
* Interior quasi-optimality examines the RMS derivative of the fitted function
  with respect to log regularization. It accepts a strict interior local
  minimum. Constant paths, the oversmoothed limit and numerical convergence
  plateaus are inadmissible; absent valleys are explicit failures. This is an
  interior version of the heuristic, not minimization over an unbounded path
  where the derivative can vanish at an endpoint.

All choices require an interior point with successful neighboring fits.
Boundary choices are recorded and the grid is extended; unresolved cases raise
structured errors. No zero penalties, default-selector fallbacks, prediction
clipping, population functions or held-out estimator observations enter the
choice. During ensemble cross-validation each path uses only that learner's
training rows; the final refit selects on its full supplied training sample.
The path is stored as `penalty_path` (also `penalty_cv` for compatibility),
with scores, residuals, seminorms, trace diagnostics, failures and the selected
point. `penalty_boundary` retains failed boundary attempts.

These rules select regularization for the empirical inverse equation. They
do not directly estimate full-function MSE and do not remove nonuniqueness.
GCV's classical fixed-operator white-noise assumptions need not hold in our
random, heteroscedastic conditional-moment problem; its nonlinear versions are
local or numerical sensitivities. L-curve and quasi-optimality are heuristics.
Candidate classes, instrument bases and the bridge link are kept fixed in the
comparison. Some individual candidates need not contain a true solution even
when the library as a whole does. These checks supply no new causal coverage
results.

Primary references: [Golub, Heath and Wahba (1979)](https://doi.org/10.1080/00401706.1979.10489751),
[Hansen and O'Leary (1993)](https://doi.org/10.1137/0914086),
[Hansen's Regularization Tools, Section 2.9](https://www2.imm.dtu.dk/~pcha/Regutools/RTv4manual.pdf),
and [Bauer and Reiss](https://arxiv.org/abs/0710.1045).
