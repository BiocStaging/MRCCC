# MRCCC 0.99.1

* The priors on the second-stage coefficients `(beta_X, beta_XZ)` and
  `beta_Z` are now scaled by a fixed diagonal computed once from the
  least-squares first stage, instead of by the current `X*`, `Z*` at every
  sweep. Because `X*` and `Z*` are functions of the first-stage parameters,
  the earlier prior entered their full conditionals, which the sampler did not
  account for. With the scale fixed the first-stage updates are exact as
  written and every update remains conjugate. The scale is returned in
  `settings$prior_scale`. `mr_ccc_gibbs()` gains `legacy_latent_prior`
  (default `FALSE`) to reproduce earlier fits exactly.

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
