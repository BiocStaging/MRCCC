#' Fit the MR-CCC model to one ligand-receptor-pathway triplet
#'
#' Fits the Bayesian Mendelian randomization model for causal cell-cell
#' communication to a single triplet consisting of ligand expression in a
#' sender cell type, receptor expression in a receiver cell type, and
#' pathway activity in the receiver cell type, using cis-eQTL genotypes as
#' instruments. The function validates and centres the inputs, runs one or
#' more Gibbs chains through [mr_ccc_gibbs()], pools the retained draws, and
#' computes the posterior inclusion probability (PIP) together with its Monte
#' Carlo standard error and convergence diagnostics.
#'
#' @param X Ligand expression in the sender cell type: a numeric vector of
#'   length \eqn{n} or an \eqn{n \times 1} matrix, one value per donor. A
#'   `SummarizedExperiment` with exactly one row (gene) is accepted directly;
#'   its first assay row is used. [se_row()] selects a gene from a larger
#'   `SummarizedExperiment`.
#' @param Z Receptor expression in the receiver cell type, same format as
#'   `X` (including a single-row `SummarizedExperiment`).
#' @param Y Pathway activity in the receiver cell type, same format as `X`
#'   (including a single-row `SummarizedExperiment`).
#' @param G Numeric matrix with \eqn{n} rows: cis-eQTL genotypes that
#'   instrument the ligand (at least one column).
#' @param H Numeric matrix with \eqn{n} rows: cis-eQTL genotypes that
#'   instrument the receptor (at least one column).
#' @param V Numeric matrix with \eqn{n} rows of shared donor covariates
#'   (for example genotype principal components, age or sex). At least one
#'   column is required; supplying `NULL` is an error. When no covariates are
#'   available, pass a single centred, non-constant column such as
#'   standardised age.
#' @param n_iter Total number of Gibbs iterations per chain.
#' @param burn_in Number of initial iterations discarded per chain. Must be
#'   smaller than `n_iter`.
#' @param thin Thinning interval: every `thin`-th post-burn-in iteration is
#'   retained.
#' @param n_chains Number of independent chains, pooled for all posterior
#'   summaries. With two or more chains the Gelman-Rubin statistic is reported
#'   and the chains start from dispersed values. Pooling `C` chains of `n_iter`
#'   iterations gives the same Monte Carlo precision as one chain of
#'   `C * n_iter`, because per-chain standard errors combine in quadrature, so
#'   raising `n_chains` is an alternative to raising `n_iter` rather than an
#'   additional cost.
#' @param init_scale Standard deviation of the normal distribution from which
#'   \eqn{\beta_X} and \eqn{\beta_{XZ}} are initialised when `n_chains > 1`.
#'   With a single chain the sampler starts deterministically at zero.
#' @param pip_threshold Threshold in \eqn{(0, 1)} that defines a discovery
#'   (`PIP > pip_threshold`). The default of 0.5 corresponds to the median
#'   probability model. Used only for the borderline warning and for the
#'   summary methods.
#' @param hyper Prior hyperparameters, as returned by
#'   [mrccc_hyperparameters()].
#' @param seed Optional integer; if supplied, the random number generator is
#'   seeded once before the first chain, which makes the whole multi-chain fit
#'   reproducible as a unit. The default `NULL` leaves the generator untouched,
#'   so a caller who prefers to manage the random number stream themselves can
#'   simply call [set.seed()] beforehand and omit this argument.
#' @param verbose Logical; if `TRUE`, each chain prints a progress line every
#'   1000 iterations.
#'
#' @details
#' **Model.** With \eqn{X^* = E[X \mid G, V]} and \eqn{Z^* = E[Z \mid H, V]}
#' the genetically predicted ligand and receptor levels, the outcome model is
#' \deqn{Y = \mu + \gamma(\beta_X X^* + \beta_{XZ} X^* Z^*) + \beta_Z Z^*
#'   + \alpha_Y^\top V + \varepsilon_Y,}
#' where \eqn{\gamma \in \{0, 1\}} is a spike-and-slab inclusion indicator
#' for the causal block \eqn{(\beta_X, \beta_{XZ})}. The indicator is
#' implemented through the prior variance of that block (a continuous
#' spike-and-slab: under \eqn{\gamma = 0} the prior variance is multiplied by
#' `nu1`, \eqn{10^{-4}} by default), so that under the spike the two
#' coefficients are concentrated tightly at zero rather than fixed at zero;
#' at this value of `nu1` the two formulations are numerically
#' indistinguishable. The posterior inclusion probability
#' \eqn{\mathrm{PIP} = P(\gamma = 1 \mid \mathrm{data})} is the
#' communication score. The receptor-modulated ligand effect is
#' \eqn{\beta_X + \beta_{XZ} z}, which changes sign at
#' \eqn{\tau = -\beta_X / \beta_{XZ}}; see [sign_reversal()].
#'
#' **Centring.** The model requires `X`, `Z` and `Y` to have mean zero.
#' They are centred internally; the standard deviations of the original
#' inputs are stored in the `scale` element so that effects can be reported
#' on a standardised scale (per standard deviation of \eqn{X} and \eqn{Z}, in
#' standard deviations of \eqn{Y}).
#'
#' **Priors.** The first-stage blocks \eqn{\pi_X, \alpha_X, \pi_Z, \alpha_Z}
#' and the covariate block \eqn{\alpha_Y} carry Zellner g-priors on their
#' observed designs \eqn{G}, \eqn{H} and \eqn{V}. The second-stage blocks are
#' scaled by the size of their regressors, but with that scale fixed:
#' \deqn{(\beta_X, \beta_{XZ}) \mid \gamma, \sigma_Y^2 \sim
#'   N_2\!\left(0,\; g\, s_\gamma\, \sigma_Y^2\, D_0^{-1}\right), \qquad
#'   \beta_Z \mid \sigma_Y^2 \sim N\!\left(0,\; g\, \sigma_Y^2 / d_Z\right),}
#' where \eqn{D_0 = \mathrm{diag}(\|\hat X^*\|^2, \|\hat X^* \circ \hat Z^*\|^2)}
#' and \eqn{d_Z = \|\hat Z^*\|^2} are computed once from the least-squares
#' first stage and held fixed, and \eqn{s_\gamma} is 1 in the slab and `nu1`
#' in the spike. Equivalently, the causal coefficients carry independent
#' normal priors on the scale of the standardised regressors. The scale is
#' fixed rather than recomputed from the current draws because \eqn{X^*} and
#' \eqn{Z^*} are functions of the first-stage parameters; a prior whose
#' covariance moved with them would enter their full conditionals, and the
#' sampler would no longer be a Gibbs sampler for the stated model. The
#' values used are returned in `settings$prior_scale`. Unless overridden in
#' `hyper`, every scale \eqn{g} is set to \eqn{\min(n, 100)}.
#'
#' **Instrument strength.** The first-stage partial F statistic of each
#' instrument block, after the covariates, is computed by ordinary least
#' squares and returned in `instruments`. When either F is below 10 a message
#' is issued. The threshold is the conventional two-stage least-squares
#' heuristic and is reported for orientation: the one-sample Bayesian model
#' fitted here does not require it, but posterior inclusion probabilities are
#' less decisive under weak instruments and longer chains may be needed.
#'
#' **Chains and pooling.** Chain \eqn{c} starts with \eqn{\gamma = 1} when
#' \eqn{c} is odd and \eqn{\gamma = 0} when \eqn{c} is even, and, when
#' `n_chains > 1`, with \eqn{(\beta_X, \beta_{XZ})} drawn from a normal
#' distribution with mean zero and standard deviation `init_scale`.
#' Retained draws are pooled across
#' chains for all posterior summaries.
#'
#' **Monte Carlo error of the PIP.** The PIP is the mean of an autocorrelated
#' binary chain, so its Monte Carlo standard error (MCSE) is estimated by
#' non-overlapping batch means within each chain and combined across chains
#' in quadrature: \eqn{\mathrm{MCSE} = \sqrt{\sum_c \mathrm{MCSE}_c^2} / C}.
#' Effective sample sizes add across chains. Batch means are never applied to
#' the concatenated chains, because a batch straddling the join between two
#' chains started at different values would inflate the variance estimate.
#' If \eqn{|\mathrm{PIP} - t| < 2\,\mathrm{MCSE}}, where \eqn{t} is
#' `pip_threshold`, the
#' discovery status of the triplet is not resolved at the current chain
#' length and a warning is issued.
#'
#' **Convergence diagnostics.** The Gelman-Rubin statistic \eqn{\hat R} is
#' reported for \eqn{\beta_X}, \eqn{\beta_{XZ}}, \eqn{\gamma} and the joint
#' log-likelihood. It compares the variance between chains with the variance
#' within them, so it requires `n_chains >= 2` and is `NA` for a single chain;
#' it is also only informative when the chains start from dispersed values,
#' which `init_scale > 0` and the alternating \eqn{\gamma} start provide.
#' Values below 1.01 are conventionally taken as consistent with convergence.
#' Because a useful \eqn{\hat R} is the reason to run several chains at all,
#' the default `n_chains` is 4.
#'
#' @return An object of class `mrccc_fit`: a list with elements
#' \describe{
#'   \item{`pip`}{Posterior inclusion probability, pooled over chains.}
#'   \item{`pip_mcse`}{Monte Carlo standard error of `pip`.}
#'   \item{`pip_ess`}{Effective sample size of the pooled \eqn{\gamma}
#'     chain.}
#'   \item{`estimates`}{Data frame with columns `parameter` (`beta_X`,
#'     `beta_XZ`, `beta_Z`), `raw_mean` (posterior mean on the centred input
#'     scale) and `std_mean` (posterior mean on the standardised scale).}
#'   \item{`draws`}{List of pooled numeric vectors `beta_X`, `beta_XZ`,
#'     `beta_Z`, `gamma`, `mu` and `loglik`, the first five on the raw
#'     (centred) scale, each of length `n_keep * n_chains`. `loglik` is the
#'     joint log-likelihood of the three equations at each retained draw and
#'     is useful for trace plots.}
#'   \item{`scale`}{List with `sd_X`, `sd_Z`, `sd_Y` (standard deviations of
#'     the original inputs), `mean_X`, `mean_Z`, `mean_Y` (the subtracted
#'     means) and `Z_range_std` (range of the centred receptor expression in
#'     standard-deviation units).}
#'   \item{`diagnostics`}{Data frame with columns `parameter` (`beta_X`,
#'     `beta_XZ`, `gamma`, `loglik`) and `rhat`.}
#'   \item{`instruments`}{Data frame with one row per exposure (`X`, `Z`)
#'     and columns `n_instruments`, `F_statistic` (first-stage partial F of
#'     the instrument block after the covariates) and `weak` (`TRUE` when
#'     `F_statistic < 10`).}
#'   \item{`settings`}{List with `n_iter`, `burn_in`, `thin`, `n_chains`,
#'     `n_keep` (retained draws per chain), `pip_threshold` and `hyper`.}
#'   \item{`call`}{The matched call.}
#' }
#'
#' @seealso [summary.mrccc_fit()], [credible_intervals()],
#'   [sign_probabilities()], [sign_reversal()], [plot_effect_curve()],
#'   [bayes_fdr()], [simulate_mrccc()], [mrccc_hyperparameters()].
#'
#' @references
#' Sarkar B, Ni Y (2026). MR-CCC: Bayesian Mendelian Randomization for
#' Causal Cell-Cell Communication. arXiv:2604.23917.
#'
#' @examples
#' sim <- simulate_mrccc(n = 150, seed = 1)
#' fit <- mr_ccc(sim$X, sim$Z, sim$Y, sim$G, sim$H, sim$V,
#'               n_iter = 600, burn_in = 100, seed = 1)
#' fit
#' fit$pip
#' fit$estimates
#'
#' @export
mr_ccc <- function(X, Z, Y, G, H, V = NULL,
                   n_iter = 20000, burn_in = 2000, thin = 1, n_chains = 4,
                   init_scale = 1, pip_threshold = 0.5,
                   hyper = mrccc_hyperparameters(), seed = NULL,
                   verbose = FALSE) {

  cl <- match.call()

  # ---- SummarizedExperiment coercion ----------------------------------------
  # A single-row SummarizedExperiment supplied as X, Z or Y is reduced to its
  # first assay row; larger objects require a gene selection through se_row().
  X <- coerce_se_vector(X, "X")
  Z <- coerce_se_vector(Z, "Z")
  Y <- coerce_se_vector(Y, "Y")

  # ---- Input validation -----------------------------------------------------
  X <- as_column(X, "X")
  Z <- as_column(Z, "Z")
  Y <- as_column(Y, "Y")
  G <- as_design(G, "G")
  H <- as_design(H, "H")
  if (is.null(V)) {
    stop("'V' must be a numeric matrix with at least one covariate column; ",
         "NULL is not allowed. When no covariates are available, supply a ",
         "single centred, non-constant column such as standardised age.",
         call. = FALSE)
  }
  V <- as_design(V, "V")

  n <- nrow(X)
  dims <- c(X = nrow(X), Z = nrow(Z), Y = nrow(Y),
            G = nrow(G), H = nrow(H), V = nrow(V))
  if (any(dims != n)) {
    stop("All inputs must have the same number of rows (donors); ",
         "row counts are: ",
         paste(names(dims), dims, sep = " = ", collapse = ", "), ".",
         call. = FALSE)
  }
  if (n < 10L) {
    stop("At least 10 donors are required; ", n, " supplied.", call. = FALSE)
  }
  inputs <- list(X = X, Z = Z, Y = Y, G = G, H = H, V = V)
  for (nm in names(inputs)) {
    if (any(!is.finite(inputs[[nm]]))) {
      stop("'", nm, "' contains NA, NaN or infinite values.", call. = FALSE)
    }
  }

  # The three design matrices are additionally checked for constant columns,
  # duplicated columns and saturation, each of which yields a fit that is
  # numerically defined but not interpretable.
  check_design(G, "G")
  check_design(H, "H")
  check_design(V, "V")

  # The combined first- and second-stage regressions need more donors than
  # coefficients. Checking the total here gives one clear message rather than
  # an implausible posterior.
  p_total <- ncol(G) + ncol(H) + ncol(V)
  if (n <= p_total + 4L) {
    stop("Too few donors for the model: ", n, " donors against ",
         p_total, " instrument and covariate columns plus four second-stage ",
         "coefficients. Supply more donors or fewer instruments.",
         call. = FALSE)
  }

  check_scalar(n_iter, "n_iter")
  check_scalar(burn_in, "burn_in")
  check_scalar(thin, "thin")
  check_scalar(n_chains, "n_chains")
  check_scalar(init_scale, "init_scale")
  check_scalar(pip_threshold, "pip_threshold")
  if (n_iter != round(n_iter) || n_iter < 1) {
    stop("'n_iter' must be a positive integer.", call. = FALSE)
  }
  if (burn_in != round(burn_in) || burn_in < 0) {
    stop("'burn_in' must be a non-negative integer.", call. = FALSE)
  }
  if (n_iter <= burn_in) {
    stop("'n_iter' (", n_iter, ") must be greater than 'burn_in' (",
         burn_in, ").", call. = FALSE)
  }
  if (thin != round(thin) || thin < 1) {
    stop("'thin' must be an integer greater than or equal to 1.",
         call. = FALSE)
  }
  if (n_chains != round(n_chains) || n_chains < 1) {
    stop("'n_chains' must be an integer greater than or equal to 1.",
         call. = FALSE)
  }
  if (init_scale < 0) {
    stop("'init_scale' must be non-negative.", call. = FALSE)
  }
  if (pip_threshold <= 0 || pip_threshold >= 1) {
    stop("'pip_threshold' must lie strictly between 0 and 1.", call. = FALSE)
  }
  if (!is(hyper, "mrccc_hyper")) {
    stop("'hyper' must be created by mrccc_hyperparameters().", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("'verbose' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(seed)) {
    check_scalar(seed, "seed")
    if (seed != round(seed)) {
      stop("'seed' must be a whole number.", call. = FALSE)
    }
  }

  # Retained draws per chain. Two are the minimum for any variance-based
  # summary, so a configuration that retains fewer is rejected here rather
  # than producing NA diagnostics without explanation.
  n_keep_expected <- (n_iter - burn_in) %/% thin
  if (n_keep_expected < 2L) {
    stop("The chosen 'n_iter' (", n_iter, "), 'burn_in' (", burn_in,
         ") and 'thin' (", thin, ") retain only ", n_keep_expected,
         " draw(s) per chain. At least two are required; lengthen the chain, ",
         "shorten the burn-in, or reduce the thinning interval.",
         call. = FALSE)
  }

  n_iter   <- as.integer(n_iter)
  burn_in  <- as.integer(burn_in)
  thin     <- as.integer(thin)
  n_chains <- as.integer(n_chains)

  # ---- Centring and scale -------------------------------------------------
  # The model requires X, Z and Y to have mean zero. The standard deviations
  # of the original inputs are retained for the standardised effect scale.
  sd_X <- sd(X[, 1L]); sd_Z <- sd(Z[, 1L]); sd_Y <- sd(Y[, 1L])
  if (any(!is.finite(c(sd_X, sd_Z, sd_Y))) || any(c(sd_X, sd_Z, sd_Y) <= 0)) {
    stop("'X', 'Z' and 'Y' must each have positive standard deviation.",
         call. = FALSE)
  }
  mean_X <- mean(X[, 1L]); mean_Z <- mean(Z[, 1L]); mean_Y <- mean(Y[, 1L])
  Xc <- X - mean_X
  Zc <- Z - mean_Z
  Yc <- Y - mean_Y
  Z_range_std <- range(Zc[, 1L] / sd_Z)

  # ---- Instrument strength ------------------------------------------------
  # First-stage partial F for each exposure: the F test of the instrument
  # block after the covariates, from ordinary least squares. This is the
  # conventional diagnostic of instrument strength. The value F > 10 is a
  # heuristic from two-stage least squares, not a requirement of the
  # one-sample Bayesian model fitted here; it is reported so that the regime
  # is visible, and a message is issued below when either F falls under it.
  instruments <- data.frame(
    exposure      = c("X", "Z"),
    n_instruments = c(ncol(G), ncol(H)),
    F_statistic   = c(first_stage_F(Xc[, 1L], G, V),
                      first_stage_F(Zc[, 1L], H, V)),
    stringsAsFactors = FALSE
  )
  instruments$weak <- is.finite(instruments$F_statistic) &
    instruments$F_statistic < 10

  if (any(instruments$weak)) {
    wk <- instruments[instruments$weak, , drop = FALSE]
    message("First-stage F below 10 for ",
            paste0(wk$exposure, " (F = ", vapply(wk$F_statistic, fmt, ""),
                   ")", collapse = " and "),
            ": the instruments are weak by the two-stage least-squares ",
            "heuristic. MR-CCC does not require F > 10, but posterior ",
            "inclusion probabilities are less decisive in this regime and ",
            "longer chains may be needed to resolve them; see 'F_statistic' ",
            "in the returned 'instruments' table.")
  }

  # ---- g-prior scale ------------------------------------------------------
  g_val <- if (is.null(hyper$g)) min(n, 100) else hyper$g

  # ---- Run the chains -----------------------------------------------------
  # The generator is seeded ONLY when the caller supplies a seed; with the
  # default seed = NULL it is left untouched and the chains draw from the
  # ambient stream. Seeding here rather than asking the caller to do it makes a
  # multi-chain fit reproducible as a unit: the chains run consecutively
  # through one stream, so a single seed fixes all of them.
  if (!is.null(seed)) set.seed(as.integer(seed))
  fits <- vector("list", n_chains)
  for (cc in seq_len(n_chains)) {
    fits[[cc]] <- mr_ccc_gibbs(
      Xc, Zc, Yc, G, H, V,
      n_iter = n_iter, burn_in = burn_in, thin = thin,
      a_sigma = hyper$a_sigma, b_sigma = hyper$b_sigma,
      a_rho = hyper$a_rho, b_rho = hyper$b_rho,
      nu1 = hyper$nu1,
      gG = g_val, gV = g_val, gH = g_val, gZ = g_val, gBeta = g_val,
      ridge = hyper$ridge,
      init_gamma = if (cc %% 2L == 0L) 0L else 1L,
      init_scale = if (n_chains > 1L) init_scale else 0,
      verbose = verbose
    )
  }
  n_keep <- as.integer(fits[[1L]]$n_keep)

  # ---- Per-chain and pooled draws -------------------------------------------
  chain_of <- function(name) lapply(fits, function(f) as.numeric(f[[name]]))
  bX_ch  <- chain_of("Beta_X_draws")
  bXZ_ch <- chain_of("Beta_XZ_draws")
  bZ_ch  <- chain_of("Beta_Z_draws")
  g_ch   <- chain_of("gamma_draws")
  mu_ch  <- chain_of("mu_draws")
  ll_ch  <- chain_of("loglik_draws")

  draws <- list(
    beta_X  = unlist(bX_ch,  use.names = FALSE),
    beta_XZ = unlist(bXZ_ch, use.names = FALSE),
    beta_Z  = unlist(bZ_ch,  use.names = FALSE),
    gamma   = unlist(g_ch,   use.names = FALSE),
    mu      = unlist(mu_ch,  use.names = FALSE),
    loglik  = unlist(ll_ch,  use.names = FALSE)
  )

  # ---- Posterior summaries ------------------------------------------------
  pip <- mean(draws$gamma)
  raw_mean <- c(beta_X  = mean(draws$beta_X),
                beta_XZ = mean(draws$beta_XZ),
                beta_Z  = mean(draws$beta_Z))
  sf <- c(beta_X  = sd_X / sd_Y,
          beta_XZ = sd_X * sd_Z / sd_Y,
          beta_Z  = sd_Z / sd_Y)
  estimates <- data.frame(
    parameter = names(raw_mean),
    raw_mean  = unname(raw_mean),
    std_mean  = unname(raw_mean * sf),
    stringsAsFactors = FALSE
  )

  # ---- Monte Carlo error of the PIP ---------------------------------------
  mc <- mcse_pooled(lapply(g_ch, mcse_batch))
  pip_mcse <- unname(mc[["mcse"]])
  pip_ess  <- unname(mc[["ess"]])

  # ---- Convergence diagnostics --------------------------------------------
  # The Gelman-Rubin statistic is the reported diagnostic. It compares the
  # variance between chains with the variance within them, so it requires at
  # least two chains and is NA for a single chain. The joint log-likelihood is
  # included alongside the three parameters because it is a scalar every
  # parameter feeds into, so drift elsewhere in the model tends to surface in
  # it; it is not a substitute for the indicator's own diagnostic, since the
  # log-likelihood is dominated by the residual variances.
  diag_pars <- c("beta_X", "beta_XZ", "gamma", "loglik")
  diag_ch   <- list(bX_ch, bXZ_ch, g_ch, ll_ch)
  diagnostics <- data.frame(
    parameter = diag_pars,
    rhat      = if (n_chains >= 2L) {
      vapply(diag_ch, gelman_rubin, numeric(1L))
    } else {
      rep(NA_real_, length(diag_pars))
    },
    stringsAsFactors = FALSE
  )

  # ---- Borderline warning -------------------------------------------------
  # The default n_iter is deliberately modest so that a first call returns
  # quickly. When the PIP lands close enough to the threshold that Monte Carlo
  # error alone could move it across, the warning names the concrete lengths
  # that resolved the published analysis rather than advising "increase
  # n_iter" without a scale. Because the MCSE falls as the reciprocal square
  # root of the number of draws, moving from 20,000 to 100,000 roughly halves
  # it and 400,000 roughly quarters it; the suggested next step therefore
  # depends on how far the PIP currently sits from the threshold.
  if (is.finite(pip_mcse) && abs(pip - pip_threshold) < 2 * pip_mcse) {
    suggestion <- if (n_iter < 100000L) {
      "Increase 'n_iter' to 100000, or to 400000 if the PIP remains borderline"
    } else if (n_iter < 400000L) {
      "Increase 'n_iter' to 400000"
    } else {
      "Increase 'n_chains', which pools additional draws at proportional cost"
    }
    warning("The posterior inclusion probability (", fmt(pip),
            ") lies within 2 * MCSE (MCSE = ", fmt(pip_mcse),
            ") of the discovery threshold (", pip_threshold,
            "): the discovery status of this triplet is not resolved at ",
            "this chain length. ", suggestion,
            ". Pooling several chains is equivalent in precision to one ",
            "chain of the same total length and additionally yields the ",
            "Gelman-Rubin statistic.",
            call. = FALSE)
  }

  structure(
    list(
      pip         = pip,
      pip_mcse    = pip_mcse,
      pip_ess     = pip_ess,
      estimates   = estimates,
      draws       = draws,
      scale       = list(sd_X = sd_X, sd_Z = sd_Z, sd_Y = sd_Y,
                         mean_X = mean_X, mean_Z = mean_Z, mean_Y = mean_Y,
                         Z_range_std = Z_range_std),
      diagnostics = diagnostics,
      instruments = instruments,
      settings    = list(n_iter = n_iter, burn_in = burn_in, thin = thin,
                         n_chains = n_chains, n_keep = n_keep,
                         pip_threshold = pip_threshold, hyper = hyper,
                         # Fixed plug-in scale of the second-stage priors.
                         # Deterministic in the data, so identical across
                         # chains; taken from the first.
                         prior_scale = fits[[1L]]$prior_scale),
      call        = cl
    ),
    class = "mrccc_fit"
  )
}
