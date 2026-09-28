#' Plot the receptor-modulated ligand effect with a credible band
#'
#' Draws the posterior mean of the receptor-modulated ligand effect
#' \eqn{\mathrm{effect}(z) = \beta_X + \beta_{XZ} z} over the observed range
#' of receptor expression, together with a pointwise credible band computed
#' from paired posterior draws.
#'
#' @param fit An object of class `mrccc_fit` returned by [mr_ccc()].
#' @param Z_observed Optional numeric vector of observed receptor expression
#'   on the original (uncentred) scale, typically the `Z` passed to
#'   [mr_ccc()]. It determines the horizontal range of the curve. If `NULL`,
#'   the range stored in the fit is used.
#' @param level Credible level in \eqn{(0, 1)} of the pointwise band.
#' @param n_draws Maximum number of paired draws used for the band. A
#'   systematic subsample of the pooled draws of this size is taken; the
#'   default of 2000 is ample for the quantiles of the band.
#' @param standardized Logical; if `TRUE` (the default) the effect is shown on
#'   the standardised scale (standard deviations of \eqn{Y} per standard
#'   deviation of \eqn{X}) against receptor expression in standard-deviation
#'   units relative to its mean. If `FALSE`, the effect is on the centred
#'   input scale against centred receptor expression.
#'
#' @details
#' The band is built from paired draws of \eqn{(\beta_X, \beta_{XZ})}, for
#' the same reason as the sign-reversal threshold in [sign_reversal()]: the
#' curve depends on two correlated parameters at once, so combining their
#' marginal intervals would misstate the uncertainty. At each grid point
#' \eqn{z}, the band is the equal-tailed interval of
#' \eqn{\beta_X^{(s)} + \beta_{XZ}^{(s)} z} across the subsampled draws. A
#' horizontal reference line marks an effect of zero; where the band crosses
#' it, the direction of the ligand effect is uncertain.
#'
#' The draws include both spike (\eqn{\gamma = 0}) and slab
#' (\eqn{\gamma = 1}) iterations, so the band reflects the full spike-and-
#' slab posterior.
#'
#' @return A `ggplot` object that can be further modified.
#'
#' @seealso [mr_ccc()], [sign_reversal()], [plot.mrccc_fit()]
#'
#' @examples
#' sim <- simulate_mrccc(n = 150, beta_X = 0.3, beta_XZ = -0.3, seed = 1)
#' fit <- mr_ccc(sim$X, sim$Z, sim$Y, sim$G, sim$H, sim$V,
#'               n_iter = 600, burn_in = 100, seed = 1)
#' p <- plot_effect_curve(fit, Z_observed = sim$Z)
#' p
#'
#' @export
plot_effect_curve <- function(fit, Z_observed = NULL, level = 0.95,
                              n_draws = 2000, standardized = TRUE) {
  check_fit(fit)
  check_scalar(level, "level")
  check_scalar(n_draws, "n_draws")
  if (level <= 0 || level >= 1) {
    stop("'level' must lie strictly between 0 and 1.", call. = FALSE)
  }
  if (n_draws < 2 || n_draws != round(n_draws)) {
    stop("'n_draws' must be an integer greater than or equal to 2.",
         call. = FALSE)
  }
  if (!is.logical(standardized) || length(standardized) != 1L ||
      is.na(standardized)) {
    stop("'standardized' must be TRUE or FALSE.", call. = FALSE)
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
  if (!standardized) z_rng <- z_rng * sd_Z

  # Paired draws on the requested scale
  sf <- if (standardized) std_factors(fit) else c(beta_X = 1, beta_XZ = 1,
                                                  beta_Z = 1)
  bx  <- fit$draws$beta_X  * sf[["beta_X"]]
  bxz <- fit$draws$beta_XZ * sf[["beta_XZ"]]
  n_tot <- length(bx)
  if (n_tot > n_draws) {
    idx <- round(seq(1, n_tot, length.out = n_draws))
    bx  <- bx[idx]
    bxz <- bxz[idx]
  }

  # Effect curve on a grid, one row per draw
  zg <- seq(z_rng[1L], z_rng[2L], length.out = 101L)
  curves <- outer(bx, rep(1, length(zg))) + outer(bxz, zg)
  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
  band <- apply(curves, 2L, quantile, probs = probs, names = FALSE)

  e <- fit$estimates
  mean_bx  <- e$raw_mean[e$parameter == "beta_X"]  * sf[["beta_X"]]
  mean_bxz <- e$raw_mean[e$parameter == "beta_XZ"] * sf[["beta_XZ"]]

  df <- data.frame(
    z      = zg,
    effect = mean_bx + mean_bxz * zg,
    lower  = band[1L, ],
    upper  = band[2L, ]
  )

  x_lab <- if (standardized) {
    expression(paste("Receptor expression ", (Z - bar(Z)) / s[Z]))
  } else {
    expression(paste("Receptor expression ", Z - bar(Z)))
  }
  y_lab <- if (standardized) {
    expression(paste("Ligand effect ",
                     hat(beta)[X]^{(s)} + hat(beta)[XZ]^{(s)} * z))
  } else {
    expression(paste("Ligand effect ", hat(beta)[X] + hat(beta)[XZ] * z))
  }

  ggplot(df, aes(x = .data$z)) +
    geom_ribbon(aes(ymin = .data$lower, ymax = .data$upper),
                fill = "grey70", alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_line(aes(y = .data$effect), colour = "black", linewidth = 0.8) +
    labs(x = x_lab, y = y_lab,
         title = "Receptor-modulated ligand effect",
         subtitle = sprintf("Posterior mean with %d%% pointwise credible band",
                            as.integer(round(100 * level)))) +
    theme_bw()
}
