#' Simulate one ligand-receptor-pathway triplet from the MR-CCC working model
#'
#' Generates donor-level data for a single triplet under the MR-CCC data
#' generating process, including cis-eQTL instruments for the ligand and the
#' receptor, observed covariates, and an unmeasured confounder that loads on
#' the ligand, the receptor and the pathway. The output can be passed
#' directly to [mr_ccc()].
#'
#' @param n Number of donors.
#' @param pG Number of sender instruments (columns of `G`).
#' @param pH Number of receiver instruments (columns of `H`).
#' @param pV Number of observed covariates (columns of `V`); at least one.
#' @param beta_X True ligand main effect on the pathway.
#' @param beta_XZ True ligand-by-receptor interaction effect.
#' @param beta_Z True receptor main effect on the pathway.
#' @param pi_val Common first-stage effect of every instrument on its
#'   target gene (instrument strength).
#' @param conf_strength Loading of the unmeasured confounder \eqn{U} on
#'   \eqn{X}, \eqn{Z} and \eqn{Y}. Zero removes confounding.
#' @param seed Optional integer; if supplied, the data are simulated under
#'   this seed and the caller's random number state is restored afterwards.
#'   The default `NULL` uses the current random number stream.
#'
#' @details
#' With \eqn{G_{ij}, H_{ij}, V_{ij}, U_i \sim N(0, 1)} independently, the
#' model is
#' \deqn{X_i = G_i^\top \pi_X + V_i^\top \alpha_X + c\,U_i + \varepsilon_{Xi},}
#' \deqn{Z_i = H_i^\top \pi_Z + V_i^\top \alpha_Z + c\,U_i + \varepsilon_{Zi},}
#' \deqn{Y_i = \beta_X X_i + \beta_Z Z_i + \beta_{XZ} X_i Z_i
#'   + V_i^\top \alpha_Y + c\,U_i + \varepsilon_{Yi},}
#' where \eqn{c} is `conf_strength`, every element of \eqn{\pi_X} and
#' \eqn{\pi_Z} equals `pi_val`, every element of \eqn{\alpha_X},
#' \eqn{\alpha_Z} and \eqn{\alpha_Y} equals 0.3, and the residuals are
#' standard normal. The confounder \eqn{U} induces a spurious association
#' between the exposures and the outcome that a naive regression absorbs
#' into the causal coefficients; the instruments allow MR-CCC to separate
#' the causal signal from it. Setting `beta_X = beta_XZ = 0` gives a null
#' triplet.
#'
#' @return A list with elements `X`, `Z`, `Y` (numeric matrices of dimension
#'   \eqn{n \times 1}), `G`, `H`, `V` (numeric matrices with \eqn{n} rows),
#'   `U` (the unmeasured confounder, for reference only) and `truth`, a list
#'   of the generating parameters `beta_X`, `beta_XZ`, `beta_Z`, `pi_val`,
#'   `conf_strength` and `tau` (the true sign-reversal threshold
#'   \eqn{-\beta_X / \beta_{XZ}} in raw units of \eqn{Z}, or `NA` when
#'   `beta_XZ` is zero).
#'
#' @seealso [mr_ccc()]
#'
#' @examples
#' sim <- simulate_mrccc(n = 200, seed = 42)
#' str(sim$truth)
#' dim(sim$G)
#'
#' # A null triplet
#' null_sim <- simulate_mrccc(n = 200, beta_X = 0, beta_XZ = 0, seed = 42)
#'
#' @export
simulate_mrccc <- function(n = 300, pG = 5, pH = 5, pV = 3,
                           beta_X = 0.3, beta_XZ = 0.3, beta_Z = 0.5,
                           pi_val = 0.5, conf_strength = 0.7, seed = NULL) {
  sizes <- list(n = n, pG = pG, pH = pH, pV = pV)
  for (nm in names(sizes)) {
    v <- sizes[[nm]]
    check_scalar(v, nm)
    if (v < 1 || v != round(v)) {
      stop("'", nm, "' must be a positive integer.", call. = FALSE)
    }
  }
  check_scalar(beta_X, "beta_X")
  check_scalar(beta_XZ, "beta_XZ")
  check_scalar(beta_Z, "beta_Z")
  check_scalar(pi_val, "pi_val")
  check_scalar(conf_strength, "conf_strength")
  # Seeded only when the caller asks for it; the default NULL leaves the
  # generator untouched.
  if (!is.null(seed)) {
    check_scalar(seed, "seed")
    if (seed != round(seed)) {
      stop("'seed' must be a whole number.", call. = FALSE)
    }
    withr::local_seed(as.integer(seed))
  }
  n <- as.integer(n); pG <- as.integer(pG)
  pH <- as.integer(pH); pV <- as.integer(pV)

  G <- matrix(rnorm(n * pG), n, pG)
  H <- matrix(rnorm(n * pH), n, pH)
  V <- matrix(rnorm(n * pV), n, pV)
  U <- rnorm(n)

  pi_X    <- rep(pi_val, pG)
  pi_Z    <- rep(pi_val, pH)
  alpha_X <- rep(0.3, pV)
  alpha_Z <- rep(0.3, pV)
  alpha_Y <- rep(0.3, pV)

  X <- as.numeric(G %*% pi_X + V %*% alpha_X + conf_strength * U + rnorm(n))
  Z <- as.numeric(H %*% pi_Z + V %*% alpha_Z + conf_strength * U + rnorm(n))
  Y <- as.numeric(beta_X * X + beta_Z * Z + beta_XZ * (X * Z) +
                    V %*% alpha_Y + conf_strength * U + rnorm(n))

  list(
    X = matrix(X, ncol = 1L),
    Z = matrix(Z, ncol = 1L),
    Y = matrix(Y, ncol = 1L),
    G = G, H = H, V = V, U = U,
    truth = list(beta_X = beta_X, beta_XZ = beta_XZ, beta_Z = beta_Z,
                 pi_val = pi_val, conf_strength = conf_strength,
                 tau = if (beta_XZ != 0) -beta_X / beta_XZ else NA_real_)
  )
}
