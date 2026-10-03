# cmbridge 0.2.0

- Added `fit_bridge_ensemble()` and `fit_adjoint_ensemble()` for inner
  cross-validated convex stacking with a common RKHS moment criterion.
- Added simplex optimization, exact blocked Gaussian scoring and scalable
  common Nystrom scoring, candidate-weight diagnostics, and full-sample refits.
- Adjoint ensembles score complete cases and accept a loading callback for
  fold-specific, common loading estimation.
- Added equation-level, optimization, fold-isolation, missing-data, and RNG tests.
