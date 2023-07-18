mcmcmcglmm v0.9.2 (2023-07-18)
=========================

### BUG FIXES

 * `make.mini.chains` now correctly handles one dimensional input variables.
 * `extract.parameters` now correctly calculates the priors for the residuals.
 * `extract.parameters` now correctly handles models with either no residuals or random terms.

mcmcmcglmmm v0.9.0 (2022-08-08)
=========================

### NEW FEATURES

  * First release! Contains the following functions: `extract.parameters`, `flat.prior`, `make.mini.chains`, `run.mini.chains`, `diagnose.mini.chains` and `combine.mini.chains`.
  * `run.mini.chains` now contains the `randomised.factors` option to create some null model results (by randomizing the levels of factors).

### MINOR IMPROVEMENTS

  * `mini.chains` structure has now been optimized for lower RAM footprint.
