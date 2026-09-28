# MRCCC 0.99.0

* Initial Bioconductor submission.
* `mr_ccc()` fits the MR-CCC model to one ligand-receptor-pathway triplet
  with one or more Gibbs chains, returning the posterior inclusion
  probability with its Monte Carlo standard error, effect estimates on the
  raw and standardised scales, pooled posterior draws, and Geweke and
  Gelman-Rubin diagnostics.
* `credible_intervals()`, `sign_probabilities()`, `sign_reversal()` and
  `plot_effect_curve()` summarise the posterior draws.
* `bayes_fdr()` and `bfdr_threshold_set()` provide Bayesian false discovery
  rate control across many triplets.
* `simulate_mrccc()` generates data from the working model with an
  unmeasured confounder.
* `mr_ccc_gibbs()` exposes the compiled RcppArmadillo sampler.
