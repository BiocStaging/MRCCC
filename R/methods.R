#' Print an MR-CCC fit
#'
#' Prints a compact description of an `mrccc_fit` object: the chain
#' configuration, the posterior inclusion probability with its Monte Carlo
#' standard error, and the standardised posterior means of the three effect
#' sizes.
#'
#' @param x An object of class `mrccc_fit` returned by [mr_ccc()].
#' @param digits Number of significant digits used for printing.
#' @param ... Further arguments; ignored.
#'
#' @return `x`, invisibly.
#'
#' @seealso [mr_ccc()], [summary.mrccc_fit()]
#'
#' @examples
#' sim <- simulate_mrccc(n = 150, seed = 1)
#' fit <- mr_ccc(sim$X, sim$Z, sim$Y, sim$G, sim$H, sim$V,
#'               n_iter = 600, burn_in = 100, seed = 1)
#' print(fit)
#'
#' @export
print.mrccc_fit <- function(x, digits = 3L, ...) {
  s <- x$settings
  cat("MR-CCC fit\n")
  cat("  Chains: ", s$n_chains, "; retained draws per chain: ", s$n_keep,
      " (n_iter = ", s$n_iter, ", burn_in = ", s$burn_in,
      ", thin = ", s$thin, ")\n", sep = "")
  cat("  Posterior inclusion probability: ", fmt(x$pip, digits),
      "  (MCSE ", fmt(x$pip_mcse, digits), ", ESS ",
      fmt(x$pip_ess, digits), ")\n", sep = "")
  if (!is.null(x$instruments)) {
    fi <- x$instruments
    iX <- fi$exposure == "X"
    iZ <- fi$exposure == "Z"
    cat("  First-stage F: X = ", fmt(fi$F_statistic[iX], digits),
        " (", fi$n_instruments[iX], " instruments), Z = ",
        fmt(fi$F_statistic[iZ], digits),
        " (", fi$n_instruments[iZ], " instruments)",
        if (any(fi$weak)) "  [weak: F < 10]" else "", "\n", sep = "")
  }
  e <- x$estimates
  cat("  Standardised posterior means: beta_X = ",
      fmt(e$std_mean[e$parameter == "beta_X"], digits),
      ", beta_XZ = ", fmt(e$std_mean[e$parameter == "beta_XZ"], digits),
      ", beta_Z = ", fmt(e$std_mean[e$parameter == "beta_Z"], digits),
      "\n", sep = "")
  invisible(x)
}

#' Summarise an MR-CCC fit
#'
#' Prints the posterior inclusion probability with its Monte Carlo standard
#' error and effective sample size, the standardised posterior means and
#' credible intervals of the three effect sizes, the convergence diagnostics,
#' and a note when the discovery status is not resolved at the current chain
#' length. The assembled table is returned invisibly.
#'
#' @param object An object of class `mrccc_fit` returned by [mr_ccc()].
#' @param level Credible level for the reported intervals.
#' @param digits Number of significant digits used for printing.
#' @param ... Further arguments; ignored.
#'
#' @return Invisibly, a data frame with one row per quantity (`beta_X`,
#'   `beta_XZ`, `beta_Z`, `pip`) and columns `raw_mean`, `std_mean`,
#'   `std_lower`, `std_upper` and `rhat`. For the `pip` row the
#'   raw and standardised means both equal the PIP, the interval columns are
#'   `NA`, and `rhat` refers to the \eqn{\gamma} chain.
#'
#' @details
#' The discovery status is considered unresolved when
#' \eqn{|\mathrm{PIP} - t| < 2\,\mathrm{MCSE}}, where \eqn{t} is the
#' discovery threshold stored in the fit; in that case a note recommends a
#' longer run or additional chains. The Gelman-Rubin column is `NA` unless the
#' fit used at least two chains.
#'
#' @seealso [mr_ccc()], [credible_intervals()]
#'
#' @examples
#' sim <- simulate_mrccc(n = 150, seed = 1)
#' fit <- mr_ccc(sim$X, sim$Z, sim$Y, sim$G, sim$H, sim$V,
#'               n_iter = 600, burn_in = 100, seed = 1)
#' tab <- summary(fit)
#' tab
#'
#' @export
summary.mrccc_fit <- function(object, level = 0.95, digits = 3L, ...) {
  check_fit(object)
  ci <- credible_intervals(object, level = level, standardized = TRUE)
  e  <- object$estimates
  d  <- object$diagnostics
  idx <- match(c("beta_X", "beta_XZ"), d$parameter)
  tab <- data.frame(
    parameter = c(e$parameter, "pip"),
    raw_mean  = c(e$raw_mean, object$pip),
    std_mean  = c(e$std_mean, object$pip),
    std_lower = c(ci$lower, NA_real_),
    std_upper = c(ci$upper, NA_real_),
    # beta_Z has no diagnostic row of its own, hence the NA in third position;
    # the fourth row is the PIP, which takes the gamma chain's R-hat.
    rhat      = c(d$rhat[idx], NA_real_,
                  d$rhat[d$parameter == "gamma"]),
    stringsAsFactors = FALSE
  )

  s <- object$settings
  cat("MR-CCC fit summary\n")
  cat("  Chains: ", s$n_chains, "; retained draws per chain: ", s$n_keep,
      "\n", sep = "")
  cat("  Posterior inclusion probability: ", fmt(object$pip, digits),
      "  (MCSE ", fmt(object$pip_mcse, digits), ", ESS ",
      fmt(object$pip_ess, digits), ")\n", sep = "")
  cat("  Discovery threshold: ", s$pip_threshold, "\n", sep = "")
  if (!is.null(object$instruments)) {
    fi <- object$instruments
    cat("  First-stage F: X = ",
        fmt(fi$F_statistic[fi$exposure == "X"], digits),
        ", Z = ", fmt(fi$F_statistic[fi$exposure == "Z"], digits),
        if (any(fi$weak)) "  (weak by the F > 10 heuristic)" else "",
        "\n", sep = "")
  }
  cat("\n")

  cat("Standardised effects (posterior mean and ",
      round(100 * level), "% credible interval):\n", sep = "")
  for (i in seq_len(nrow(e))) {
    cat("  ", format(e$parameter[i], width = 8), " ",
        fmt(tab$std_mean[i], digits), "  [",
        fmt(tab$std_lower[i], digits), ", ",
        fmt(tab$std_upper[i], digits), "]\n", sep = "")
  }

  cat("\nGelman-Rubin R-hat (NA unless n_chains >= 2):\n")
  for (i in seq_len(nrow(d))) {
    cat("  ", format(d$parameter[i], width = 8), " R-hat = ",
        fmt(d$rhat[i], digits), "\n", sep = "")
  }

  if (is.finite(object$pip_mcse) &&
      abs(object$pip - s$pip_threshold) < 2 * object$pip_mcse) {
    cat("\nNote: the PIP lies within two Monte Carlo standard errors of the ",
        "discovery threshold; the discovery status is not resolved at this ",
        "chain length. Increase n_iter (100000, or 400000 if it remains ",
        "borderline) or n_chains.\n", sep = "")
  }
  invisible(tab)
}

#' Plot the receptor-modulated ligand effect of an MR-CCC fit
#'
#' Method for [plot()] that dispatches to [plot_effect_curve()].
#'
#' @param x An object of class `mrccc_fit` returned by [mr_ccc()].
#' @param Z_observed Optional numeric vector of observed receptor expression
#'   (on the original scale) used to set the plotting range. If `NULL`, the
#'   range stored in the fit is used.
#' @param ... Further arguments passed to [plot_effect_curve()].
#'
#' @return A `ggplot` object.
#'
#' @seealso [plot_effect_curve()]
#'
#' @examples
#' sim <- simulate_mrccc(n = 150, seed = 1)
#' fit <- mr_ccc(sim$X, sim$Z, sim$Y, sim$G, sim$H, sim$V,
#'               n_iter = 600, burn_in = 100, seed = 1)
#' plot(fit)
#'
#' @export
plot.mrccc_fit <- function(x, Z_observed = NULL, ...) {
  plot_effect_curve(x, Z_observed = Z_observed, ...)
}
