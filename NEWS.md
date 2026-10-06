# cmbridge 0.3.0.9021

- Allow Landweber candidates to cross-validate their positive conditioning-weight ridge on the same training-only person splits and raw U-statistic criterion as other candidates. Preserve boundary extension, failed-trial records, target classes and existing solver checks; retain the selected parameter name and full tuning history.

# cmbridge 0.3.0.9020

- Preserve batched penalty-path evaluation and warm starts. Fall back to individual trials only when a batch raises an error; retain failed-trial diagnostics and the interior minimum rule.

# cmbridge 0.2.0

- Added `fit_bridge_ensemble()` and `fit_adjoint_ensemble()` for inner
  cross-validated convex stacking with a common RKHS moment criterion.
- Added simplex optimization, exact blocked Gaussian scoring and scalable
  common Nystrom scoring, candidate-weight diagnostics, and full-sample refits.
- Adjoint ensembles score complete cases and accept a loading callback for
  fold-specific, common loading estimation.
- Added equation-level, optimization, fold-isolation, missing-data, and RNG tests.
