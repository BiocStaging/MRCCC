#' MRCCC: Bayesian Mendelian randomization for causal cell-cell communication
#'
#' MRCCC implements MR-CCC, a Bayesian Mendelian randomization framework for
#' inferring causal ligand-receptor signalling between cell types from
#' population-scale single-cell expression and genotype data. For each
#' ligand-receptor-pathway triplet, cis-eQTL genotypes instrument ligand
#' expression in the sender cell type and receptor expression in the receiver
#' cell type, and a working model with a ligand-by-receptor interaction is
#' fitted to downstream pathway activity by a blocked Gibbs sampler. A
#' spike-and-slab prior on the causal block yields a posterior inclusion
#' probability that serves as the communication score.
#'
#' @section Main functions:
#' \describe{
#'   \item{[mr_ccc()]}{Fit the model to one triplet and obtain the posterior
#'     inclusion probability, effect estimates, Monte Carlo error and
#'     convergence diagnostics.}
#'   \item{[credible_intervals()]}{Credible intervals for the three effect
#'     sizes.}
#'   \item{[sign_probabilities()]}{Marginal and conditional posterior sign
#'     probabilities of the interaction effect.}
#'   \item{[sign_reversal()]}{Posterior of the receptor level at which the
#'     ligand effect reverses sign.}
#'   \item{[plot_effect_curve()]}{Effect curve with a pointwise credible
#'     band.}
#'   \item{[bayes_fdr()], [bfdr_threshold_set()]}{Bayesian false discovery
#'     rate control across many triplets.}
#'   \item{[simulate_mrccc()]}{Simulate a triplet from the working model.}
#'   \item{[mr_ccc_gibbs()]}{The compiled Gibbs sampler.}
#' }
#'
#' @references
#' Sarkar B, Ni Y (2026). MR-CCC: Bayesian Mendelian Randomization for
#' Causal Cell-Cell Communication. arXiv:2604.23917.
#'
#' @useDynLib MRCCC, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom SummarizedExperiment assay assayNames
#' @importFrom methods is
#' @importFrom stats var sd quantile rnorm lm.fit
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_hline geom_line labs
#'   theme_bw .data
#' @keywords internal
"_PACKAGE"
