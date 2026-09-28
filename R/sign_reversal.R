#' Receptor level at which the ligand effect reverses sign
#'
#' Computes the posterior distribution of the sign-reversal threshold
#' \eqn{\tau = -\beta_X / \beta_{XZ}}, the receptor expression level at which
#' the receptor-modulated ligand effect \eqn{\beta_X + \beta_{XZ} z} changes
#' sign, from paired posterior draws of an MR-CCC fit.
#'
#' @param fit An object of class `mrccc_fit` returned by [mr_ccc()].
#' @param Z_observed Optional numeric vector of observed receptor expression
#'   on the original (uncentred) scale, typically the `Z` passed to
#'   [mr_ccc()]. If supplied, the posterior probability that \eqn{\tau} falls
#'   inside the observed range is also returned. If `NULL`, the range stored
#'   in the fit is used.
#' @param level Credible level in \eqn{(0, 1)} for the reported interval.
#'
#' @details
#' \eqn{\tau} is a ratio of two parameters, so its posterior can only be
#' obtained from paired draws of \eqn{(\beta_X, \beta_{XZ})}; combining the
#' two marginal summaries would misstate the uncertainty because the two
#' coefficients are strongly correlated a posteriori. Draws in which
#' \eqn{|\beta_{XZ}| < 10^{-8}} are dropped, since \eqn{\tau} diverges there.
#'
#' **Scale.** [mr_ccc()] centres \eqn{Z} before fitting, so on the raw scale
#' \eqn{\tau_{\mathrm{raw}} = -\beta_X / \beta_{XZ}} is expressed in the
#' original units of \eqn{Z} *relative to the mean receptor expression*.
#' The function reports \eqn{\tau = \tau_{\mathrm{raw}} / \mathrm{sd}(Z)},
#' in standard-deviation units of \eqn{Z} relative to its mean, which is the
#' scale on which the standardised coefficients are defined:
#' \eqn{-\beta_X^{(s)} / \beta_{XZ}^{(s)} = -(\beta_X / \beta_{XZ}) /
#' \mathrm{sd}(Z)}. A value of \eqn{\tau = 0} therefore corresponds to a
#' donor with average receptor expression, and \eqn{\tau = 1} to a donor one
#' standard deviation above the average. The observed range is expressed on
#' the same scale, \eqn{(Z - \bar Z) / \mathrm{sd}(Z)}.
#'
#' Whether the reversal is realised in the population is determined by
#' `p_in_range`, the posterior probability that \eqn{\tau} lies inside the
#' observed standardised receptor range, rather than by \eqn{\tau} itself:
#' a threshold far outside the observed range describes a direction change
#' that no donor exhibits. The posterior of a ratio is heavy-tailed, so the
#' median and quantiles are reported rather than a mean.
#'
#' @return A list with elements
#' \describe{
#'   \item{`tau_median`}{Posterior median of \eqn{\tau} in standard-deviation
#'     units of \eqn{Z}.}
#'   \item{`tau_lower`, `tau_upper`}{Equal-tailed credible bounds at `level`.}
#'   \item{`level`}{The credible level.}
#'   \item{`n_draws`}{Number of paired draws used.}
#'   \item{`Z_range_std`}{The observed receptor range in standard-deviation
#'     units.}
#'   \item{`p_in_range`}{Posterior probability that \eqn{\tau} lies inside
#'     `Z_range_std`.}
#'   \item{`tau_draws`}{The paired draws of \eqn{\tau} (standardised scale).}
#' }
#' All numeric summaries are `NA` when no usable paired draw is available.
#'
#' @seealso [mr_ccc()], [plot_effect_curve()], [sign_probabilities()]
#'
#' @examples
#' sim <- simulate_mrccc(n = 150, beta_X = 0.3, beta_XZ = -0.3, seed = 1)
#' fit <- mr_ccc(sim$X, sim$Z, sim$Y, sim$G, sim$H, sim$V,
#'               n_iter = 600, burn_in = 100, seed = 1)
#' sr <- sign_reversal(fit, Z_observed = sim$Z)
#' c(sr$tau_median, sr$tau_lower, sr$tau_upper, sr$p_in_range)
#'
#' @export
sign_reversal <- function(fit, Z_observed = NULL, level = 0.95) {
  check_fit(fit)
  check_scalar(level, "level")
  if (level <= 0 || level >= 1) {
    stop("'level' must lie strictly between 0 and 1.", call. = FALSE)
  }
  sd_Z <- fit$scale$sd_Z
  if (is.null(Z_observed)) {
    z_rng <- fit$scale$Z_range_std
  } else {
    if (!is.numeric(Z_observed) || any(!is.finite(Z_observed))) {
      stop("'Z_observed' must be a numeric vector without missing values.",
           call. = FALSE)
    }
    z_rng <- range((as.numeric(Z_observed) - fit$scale$mean_Z) / sd_Z)
  }

  bx  <- fit$draws$beta_X
  bxz <- fit$draws$beta_XZ
  ok  <- is.finite(bx) & is.finite(bxz) & abs(bxz) >= 1e-8
  n_ok <- sum(ok)

  if (n_ok == 0L) {
    return(list(tau_median = NA_real_, tau_lower = NA_real_,
                tau_upper = NA_real_, level = level, n_draws = 0L,
                Z_range_std = z_rng, p_in_range = NA_real_,
                tau_draws = numeric(0)))
  }

  # tau on the raw (centred) scale is in original units of Z; dividing by
  # sd(Z) expresses it in standard-deviation units relative to the mean.
  tau_raw <- -bx[ok] / bxz[ok]
  tau     <- tau_raw / sd_Z

  probs <- c((1 - level) / 2, 0.5, 1 - (1 - level) / 2)
  q <- quantile(tau, probs = probs, names = FALSE)
  list(
    tau_median  = q[2L],
    tau_lower   = q[1L],
    tau_upper   = q[3L],
    level       = level,
    n_draws     = as.integer(n_ok),
    Z_range_std = z_rng,
    p_in_range  = mean(tau >= z_rng[1L] & tau <= z_rng[2L]),
    tau_draws   = tau
  )
}
