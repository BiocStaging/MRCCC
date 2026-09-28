# MRCCC

Bayesian Mendelian randomization for causal cell-cell communication.

MRCCC implements MR-CCC, a Bayesian Mendelian randomization framework for
inferring causal ligand-receptor signalling between cell types from
population-scale single-cell expression and genotype data. For each
ligand-receptor-pathway triplet, cis-eQTL genotypes instrument ligand
expression in the sender cell type and receptor expression in the receiver
cell type, and a working model with a ligand-by-receptor interaction is fitted
to downstream pathway activity by a blocked Gibbs sampler written in
RcppArmadillo. A spike-and-slab prior on the causal block yields a posterior
inclusion probability (PIP) that serves as the communication score; the
retained draws support credible intervals, posterior sign probabilities, the
receptor level at which the ligand effect reverses sign, Monte Carlo error and
convergence diagnostics, and Bayesian false discovery rate control across
many triplets.

## Installation

Once the package is available from Bioconductor:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("MRCCC")
```

The development version can be installed from GitHub:

```r
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("bitansa/MRCCC")
```

The sampler is compiled from C++ through RcppArmadillo, so a working C++
toolchain is required (Rtools on Windows, Xcode command line tools on macOS,
`build-essential` and `gfortran` on Debian-based Linux).

## Quick start

```r
library(MRCCC)

sim <- simulate_mrccc(n = 300, beta_X = 0.3, beta_XZ = 0.3, seed = 1)
fit <- mr_ccc(sim$X, sim$Z, sim$Y, sim$G, sim$H, sim$V,
              n_iter = 20000, burn_in = 2000, seed = 1)
summary(fit)
credible_intervals(fit)
sign_probabilities(fit)
sign_reversal(fit, Z_observed = sim$Z)[c("tau_median", "p_in_range")]
plot_effect_curve(fit, Z_observed = sim$Z)
bayes_fdr(c(0.95, 0.80, 0.55, 0.20))
```

Donor-level (pseudo-bulk) expression stored in a `SummarizedExperiment`
(genes as rows, donors as columns) can be passed to `mr_ccc()` directly when
it contains a single gene, or after selecting a gene with `se_row()`.

## Chain length and convergence

`mr_ccc()` runs `n_chains` chains (default 4) from dispersed starting values
and pools the retained draws. Pooling is not more expensive than lengthening a
single chain: per-chain Monte Carlo standard errors combine in quadrature, so
`C` chains of `n_iter` iterations carry the same standard error as one chain of
`C * n_iter`. What several chains add is the Gelman-Rubin statistic R-hat,
reported in `fit$diagnostics` for both causal coefficients, the inclusion
indicator and the joint log-likelihood, and undefined for a single chain.

`mr_ccc()` also reports the first-stage partial F statistic of each
instrument block in `fit$instruments` and issues a message when either falls
below 10. That threshold is the two-stage least-squares heuristic; the
one-sample Bayesian model does not require it, but posterior inclusion
probabilities are less decisive under weak instruments and longer chains may
be needed.

The posterior inclusion probability is the mean of a strongly autocorrelated
binary chain, so its Monte Carlo standard error matters whenever the PIP is
near the discovery threshold. `mr_ccc()` reports that standard error as
`fit$pip_mcse` and warns when the PIP lies within two of them of the
threshold, since the discovery status is then unresolved at the current chain
length. The default `n_iter = 20000` is chosen so that a first call returns
quickly; the analysis in the paper uses four chains of 100,000.

## Reproducing the paper

The scripts that reproduce the simulation studies and the real-data analysis
of the paper are maintained separately at
<https://github.com/bitansa/MR-CCC>. The processed data used by those
scripts are archived on Zenodo under DOI
[10.5281/zenodo.19675075](https://doi.org/10.5281/zenodo.19675075).

## Citation

Sarkar B, Ni Y (2026). MR-CCC: Bayesian Mendelian Randomization for Causal
Cell-Cell Communication. arXiv:2604.23917.
<https://arxiv.org/abs/2604.23917>

```r
citation("MRCCC")
```

## License

MIT. See the `LICENSE` file.
