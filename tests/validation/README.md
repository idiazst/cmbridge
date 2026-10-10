# Package recovery checks

The fast unit and limited-simulation checks are in `tests/testthat`.
Run them with `testthat::test_local(stop_on_failure = TRUE)`.

The scripts here are the existing bridge, adjoint, ensemble-selection and
inverse-expit recovery checks, moved from cmbridge-tests without changing their
data-generating distributions or pass criteria. They also write y=x figures.
These larger checks run in the package's GitHub workflow and can be run locally:

```sh
R CMD INSTALL --install-tests .
Rscript tests/validation/validate_bridge.R validation-artifacts/bridge
Rscript tests/validation/validate_adjoint.R validation-artifacts/adjoint
Rscript tests/validation/validate_ensemble.R validation-artifacts/ensemble
Rscript tests/validation/validate_linked_bridge.R validation-artifacts/linked-bridge
Rscript tests/validation/plot_recovery_truth.R validation-artifacts
Rscript tests/validation/plot_ensemble_truth.R validation-artifacts/ensemble
```

Longitudinal and imputation simulation studies belong to lmtp_selfcensor_sim.
Generated plots and results are CI artifacts, not package source files.
