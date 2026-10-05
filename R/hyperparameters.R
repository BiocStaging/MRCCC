#' Hyperparameters of the MR-CCC prior
#'
#' Constructs and validates the list of prior hyperparameters used by
#' [mr_ccc()]. The defaults are those used throughout the MR-CCC analyses and
#' are weakly informative for inputs on a scale of roughly unit variance:
#' with several hundred donors the likelihood then dominates the
#' inverse-gamma and beta priors. The inverse-gamma priors are not invariant
#' to the units of `X`, `Z` and `Y`, so inputs on very different scales
#' should be standardised before fitting.
#'
#' @param a_sigma,b_sigma Shape and scale of the inverse-gamma prior on the
#'   three residual variances \eqn{\sigma_X^2, \sigma_Z^2, \sigma_Y^2}. Both
#'   must be positive.
#' @param a_rho,b_rho Shape parameters of the beta prior on the prior
#'   inclusion probability \eqn{\rho}. Both must be positive.
#' @param nu1 Spike variance multiplier in \eqn{(0, 1)}. The spike component
#'   of the prior on \eqn{(\beta_X, \beta_{XZ})} has the slab covariance
#'   scaled by `nu1`, so small values enforce a near-zero spike.
#' @param g Common g-prior scale applied to every coefficient block
#'   (first-stage sender effects, first-stage receiver effects, covariate
#'   effects, receptor main effect, and the causal block). If `NULL` (the
#'   default), [mr_ccc()] sets \eqn{g = \min(n, 100)}, where \eqn{n} is the
#'   number of donors. A positive number overrides this rule.
#' @param ridge Non-negative diagonal ridge added before every matrix
#'   inversion for numerical stability.
#'
#' @return A named list of class `mrccc_hyper` with elements `a_sigma`,
#'   `b_sigma`, `a_rho`, `b_rho`, `nu1`, `g` and `ridge`.
#'
#' @details
#' Each first-stage block \eqn{\theta} with observed design \eqn{D} (the
#' genotype matrices \eqn{G}, \eqn{H} or the covariates \eqn{V}) has the
#' Zellner g-prior \eqn{\theta \sim N(0, g\,\sigma^2 (D^\top D)^{-1})}. The
#' second-stage blocks \eqn{(\beta_X, \beta_{XZ})} and \eqn{\beta_Z} use the
#' same scale \eqn{g} but against a fixed diagonal computed once from the
#' least-squares first stage, because their regressors \eqn{X^*} and
#' \eqn{Z^*} are themselves parameters of the model; see the *Priors* section
#' of [mr_ccc()]. The default \eqn{g = \min(n, 100)} follows the
#' unit-information convention for small samples and caps the prior
#' dispersion for large ones. The spike-and-slab prior on the causal block
#' uses scale \eqn{g} in the slab and \eqn{\nu_1 g} in the spike.
#'
#' @seealso [mr_ccc()]
#'
#' @examples
#' hyper <- mrccc_hyperparameters()
#' str(hyper)
#'
#' # A more diffuse slab and a fixed g-prior scale
#' mrccc_hyperparameters(nu1 = 1e-3, g = 50)
#'
#' @export
mrccc_hyperparameters <- function(a_sigma = 3, b_sigma = 2,
                                  a_rho = 3, b_rho = 1,
                                  nu1 = 1e-4, g = NULL, ridge = 1e-8) {
  check_scalar(a_sigma, "a_sigma")
  check_scalar(b_sigma, "b_sigma")
  check_scalar(a_rho,   "a_rho")
  check_scalar(b_rho,   "b_rho")
  check_scalar(nu1,     "nu1")
  check_scalar(ridge,   "ridge")
  if (a_sigma <= 0 || b_sigma <= 0) {
    stop("'a_sigma' and 'b_sigma' must be positive.", call. = FALSE)
  }
  if (a_rho <= 0 || b_rho <= 0) {
    stop("'a_rho' and 'b_rho' must be positive.", call. = FALSE)
  }
  if (nu1 <= 0 || nu1 >= 1) {
    stop("'nu1' must lie strictly between 0 and 1.", call. = FALSE)
  }
  if (ridge < 0) {
    stop("'ridge' must be non-negative.", call. = FALSE)
  }
  if (!is.null(g)) {
    check_scalar(g, "g")
    if (g <= 0) stop("'g' must be positive or NULL.", call. = FALSE)
  }
  structure(
    list(a_sigma = as.numeric(a_sigma), b_sigma = as.numeric(b_sigma),
         a_rho = as.numeric(a_rho), b_rho = as.numeric(b_rho),
         nu1 = as.numeric(nu1), g = if (is.null(g)) NULL else as.numeric(g),
         ridge = as.numeric(ridge)),
    class = "mrccc_hyper"
  )
}
