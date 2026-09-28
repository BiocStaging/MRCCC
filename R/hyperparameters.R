#' Hyperparameters of the MR-CCC prior
#'
#' Constructs and validates the list of prior hyperparameters used by
#' [mr_ccc()]. The defaults are those used throughout the MR-CCC analyses and
#' are weakly informative for expression data on a log scale.
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
#' Under the g-prior, each coefficient block \eqn{\theta} with design
#' \eqn{D} has prior \eqn{\theta \sim N(0, g\,\sigma^2 (D^\top D)^{-1})}. The
#' default \eqn{g = \min(n, 100)} follows the unit-information convention for
#' small samples and caps the prior dispersion for large ones. The spike-and-
#' slab prior on the causal block uses the same form with scale \eqn{g} in
#' the slab and \eqn{\nu_1 g} in the spike.
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
