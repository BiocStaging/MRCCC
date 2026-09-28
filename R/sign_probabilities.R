#' Posterior sign probabilities of the interaction effect
#'
#' Computes the posterior probability that the ligand-by-receptor
#' interaction coefficient \eqn{\beta_{XZ}} is positive and that it is
#' negative, both marginally and conditionally on the causal block being
#' included (\eqn{\gamma = 1}).
#'
#' @param fit An object of class `mrccc_fit` returned by [mr_ccc()].
#'
#' @details
#' The prior on \eqn{(\beta_X, \beta_{XZ})} is a spike-and-slab mixture, so
#' the marginal posterior of \eqn{\beta_{XZ}} is a mixture of two
#' components: draws taken while \eqn{\gamma = 0}, which are concentrated
#' near zero and carry no directional information, and draws taken while
#' \eqn{\gamma = 1}, which come from the slab. The two summaries therefore
#' answer different questions and both are informative.
#'
#' * The **marginal** probabilities `p_pos_marginal` and `p_neg_marginal`
#'   use every retained draw. They combine the evidence that there is an
#'   effect at all with the evidence about its direction. By the law of
#'   total probability,
#'   \eqn{P(\beta_{XZ} > 0) = \mathrm{PIP} \cdot P(\beta_{XZ} > 0 \mid \gamma = 1)
#'   + (1 - \mathrm{PIP}) \cdot P(\beta_{XZ} > 0 \mid \gamma = 0)}, and the
#'   last factor is close to one half because the spike is symmetric about
#'   zero.
#' * The **conditional** probabilities `p_pos_conditional` and
#'   `p_neg_conditional` use only the draws with \eqn{\gamma = 1}. They
#'   describe the direction of the interaction given that the ligand acts on
#'   the pathway, and sum to one (up to draws exactly at zero).
#'
#' A triplet with a high PIP and a conditional probability near one has a
#' well-determined direction; a triplet with a modest PIP can still have a
#' decisive conditional direction, in which case the marginal probability is
#' limited mainly by the inclusion uncertainty. When no retained draw has
#' \eqn{\gamma = 1}, the conditional columns are `NA`.
#'
#' @return A one-row data frame with columns `pip`, `p_pos_marginal`,
#'   `p_neg_marginal`, `p_pos_conditional`, `p_neg_conditional` and
#'   `n_conditional` (the number of pooled draws with \eqn{\gamma = 1}).
#'
#' @seealso [mr_ccc()], [credible_intervals()], [sign_reversal()]
#'
#' @examples
#' sim <- simulate_mrccc(n = 150, seed = 1)
#' fit <- mr_ccc(sim$X, sim$Z, sim$Y, sim$G, sim$H, sim$V,
#'               n_iter = 600, burn_in = 100, seed = 1)
#' sign_probabilities(fit)
#'
#' @export
sign_probabilities <- function(fit) {
  check_fit(fit)
  bxz <- fit$draws$beta_XZ
  g   <- fit$draws$gamma
  inc <- g == 1
  n_inc <- sum(inc)
  data.frame(
    pip               = fit$pip,
    p_pos_marginal    = mean(bxz > 0),
    p_neg_marginal    = mean(bxz < 0),
    p_pos_conditional = if (n_inc > 0L) mean(bxz[inc] > 0) else NA_real_,
    p_neg_conditional = if (n_inc > 0L) mean(bxz[inc] < 0) else NA_real_,
    n_conditional     = as.integer(n_inc),
    stringsAsFactors  = FALSE
  )
}
