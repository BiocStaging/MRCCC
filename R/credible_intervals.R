#' Credible intervals for the MR-CCC effect sizes
#'
#' Computes equal-tailed posterior credible intervals for \eqn{\beta_X},
#' \eqn{\beta_{XZ}} and \eqn{\beta_Z} from the pooled retained draws of an
#' MR-CCC fit.
#'
#' @param fit An object of class `mrccc_fit` returned by [mr_ccc()].
#' @param level Credible level in \eqn{(0, 1)}; the default 0.95 gives the
#'   2.5 and 97.5 percent posterior quantiles.
#' @param standardized Logical; if `TRUE` (the default) the draws are
#'   rescaled to the standardised scale (per standard deviation of \eqn{X}
#'   and \eqn{Z}, in standard deviations of \eqn{Y}) before the quantiles are
#'   taken. If `FALSE`, the intervals are on the centred input scale.
#'
#' @details
#' The intervals are marginal posterior quantiles of each coefficient. Because
#' the prior on \eqn{(\beta_X, \beta_{XZ})} is a spike-and-slab mixture, the
#' marginal posterior of each of these two coefficients is itself a mixture of
#' a component concentrated near zero (draws with \eqn{\gamma = 0}) and a
#' slab component (draws with \eqn{\gamma = 1}). An interval that covers zero
#' therefore reflects the posterior weight of the spike as much as the
#' uncertainty within the slab; see [sign_probabilities()] for a
#' decomposition.
#'
#' @return A data frame with columns `parameter` (`beta_X`, `beta_XZ`,
#'   `beta_Z`), `lower`, `median` and `upper`.
#'
#' @seealso [mr_ccc()], [sign_probabilities()]
#'
#' @examples
#' sim <- simulate_mrccc(n = 150, seed = 1)
#' fit <- mr_ccc(sim$X, sim$Z, sim$Y, sim$G, sim$H, sim$V,
#'               n_iter = 600, burn_in = 100, seed = 1)
#' credible_intervals(fit)
#' credible_intervals(fit, level = 0.9, standardized = FALSE)
#'
#' @export
credible_intervals <- function(fit, level = 0.95, standardized = TRUE) {
  check_fit(fit)
  check_scalar(level, "level")
  if (level <= 0 || level >= 1) {
    stop("'level' must lie strictly between 0 and 1.", call. = FALSE)
  }
  if (!is.logical(standardized) || length(standardized) != 1L ||
      is.na(standardized)) {
    stop("'standardized' must be TRUE or FALSE.", call. = FALSE)
  }
  pars <- c("beta_X", "beta_XZ", "beta_Z")
  sf <- if (standardized) std_factors(fit) else c(beta_X = 1, beta_XZ = 1,
                                                  beta_Z = 1)
  probs <- c((1 - level) / 2, 0.5, 1 - (1 - level) / 2)
  q <- vapply(pars, function(p) {
    quantile(fit$draws[[p]] * sf[[p]], probs = probs, names = FALSE)
  }, numeric(3))
  data.frame(
    parameter = pars,
    lower  = unname(q[1L, ]),
    median = unname(q[2L, ]),
    upper  = unname(q[3L, ]),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}
