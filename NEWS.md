# MRCCC 0.99.2

* `mr_ccc()` now centres every column of the instrument matrices `G`, `H`
  and the covariate matrix `V`, as it already did for `X`, `Z` and `Y`.
  Neither first stage of the model has an intercept, so uncentred designs
  (for example genotype dosages coded 0/1/2) forced each first-stage fit
  through the origin and shifted the projections `X*` and `Z*` and the fixed
  prior scale. The documentation of `mr_ccc_gibbs()`, which does not centre
  its inputs, now states that all six inputs must be centred.
* `mr_ccc()` stops with an explicit message when `cbind(G, V)` or
  `cbind(H, V)` is rank deficient after centring.
* `mr_ccc()` warns when any Gelman-Rubin R-hat exceeds 1.01, and issues a
  message when only one chain is run or when the Monte Carlo standard error
  of the PIP cannot be estimated.
* `mr_ccc()` and `simulate_mrccc()` restore the caller's random number state
  after running under a supplied `seed`. The fitted object records the
  g-prior scale used in `settings$g`.
* `mr_ccc_gibbs()` validates its arguments, reports numerically singular
  matrices with an informative error, and can be interrupted from R.

# MRCCC 0.99.1

* The priors on the second-stage coefficients `(beta_X, beta_XZ)` and
  `beta_Z` are now scaled by a fixed diagonal computed once from the
  least-squares first stage, instead of by the current `X*`, `Z*` at every
  sweep. Because `X*` and `Z*` are functions of the first-stage parameters,
  the earlier prior entered their full conditionals, which the sampler did not
  account for. With the scale fixed the first-stage updates are exact as
  written and every update remains conjugate. The scale is returned in
  `settings$prior_scale`.

# MRCCC 0.99.0

* Initial Bioconductor submission.
* `mr_ccc()` fits the MR-CCC model to one ligand-receptor-pathway triplet
  with one or more Gibbs chains, returning the posterior inclusion
  probability with its Monte Carlo standard error, effect estimates on the
  raw and standardised scales, pooled posterior draws, the joint
  log-likelihood of every retained draw, and the Gelman-Rubin statistic.
* `mr_ccc()` also reports first-stage instrument strength: the partial F of
  each instrument block after the covariates, with a message when either
  block is weak by the conventional two-stage least squares threshold.
* `credible_intervals()`, `sign_probabilities()`, `sign_reversal()` and
  `plot_effect_curve()` summarise the posterior draws.
* `bayes_fdr()` and `bfdr_threshold_set()` provide Bayesian false discovery
  rate control across many triplets.
* `simulate_mrccc()` generates data from the working model with an
  unmeasured confounder.
* `mr_ccc_gibbs()` exposes the compiled RcppArmadillo sampler.
