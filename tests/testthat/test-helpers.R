test_that("bayes_fdr reproduces the nested-set definition in input order", {
  expect_equal(bayes_fdr(c(0.9, 0.8, 0.5)), c(0.1, 0.15, 0.2666667),
               tolerance = 1e-6)
  # Unsorted input: the values are returned in input order.
  expect_equal(bayes_fdr(c(0.5, 0.9, 0.8)), c(0.2666667, 0.1, 0.15),
               tolerance = 1e-6)
  expect_length(bayes_fdr(numeric(0)), 0L)
  expect_error(bayes_fdr(c(0.5, NA)), "NA")
  expect_error(bayes_fdr(c(0.5, 1.2)), "\\[0, 1\\]")
})

test_that("bfdr_threshold_set selects the largest qualifying nested set", {
  pip <- c(0.9, 0.8, 0.5, 0.2)
  # bFDR of the nested sets: 0.1, 0.15, 0.2667, 0.4
  expect_equal(bfdr_threshold_set(pip, alpha = 0.2),
               c(TRUE, TRUE, FALSE, FALSE))
  expect_equal(bfdr_threshold_set(pip, alpha = 0.3),
               c(TRUE, TRUE, TRUE, FALSE))
  expect_equal(bfdr_threshold_set(pip, alpha = 0.05),
               c(FALSE, FALSE, FALSE, FALSE))
  # Input order is preserved.
  expect_equal(bfdr_threshold_set(c(0.2, 0.9, 0.5, 0.8), alpha = 0.2),
               c(FALSE, TRUE, FALSE, TRUE))
  expect_error(bfdr_threshold_set(pip, alpha = 1), "'alpha' must lie")
})

sim <- simulate_mrccc(n = 150, beta_X = 0.3, beta_XZ = -0.3, seed = 1)
fit <- suppressWarnings(
  mr_ccc(sim$X, sim$Z, sim$Y, sim$G, sim$H, sim$V,
         n_iter = 1200, burn_in = 200, seed = 1)
)

test_that("sign_probabilities returns valid probabilities", {
  sp <- sign_probabilities(fit)
  expect_s3_class(sp, "data.frame")
  expect_equal(nrow(sp), 1L)
  expect_true(sp$p_pos_marginal >= 0 && sp$p_pos_marginal <= 1)
  expect_true(sp$p_neg_marginal >= 0 && sp$p_neg_marginal <= 1)
  expect_true(sp$p_pos_marginal + sp$p_neg_marginal <= 1 + 1e-12)
  if (sp$n_conditional > 0) {
    expect_true(sp$p_pos_conditional >= 0 && sp$p_pos_conditional <= 1)
    expect_true(sp$p_neg_conditional >= 0 && sp$p_neg_conditional <= 1)
    expect_true(sp$p_pos_conditional + sp$p_neg_conditional <= 1 + 1e-12)
  }
  expect_equal(sp$pip, fit$pip)
})

test_that("credible_intervals are ordered and respect the scale", {
  ci <- credible_intervals(fit)
  expect_equal(ci$parameter, c("beta_X", "beta_XZ", "beta_Z"))
  expect_true(all(ci$lower <= ci$median))
  expect_true(all(ci$median <= ci$upper))
  ci_raw <- credible_intervals(fit, standardized = FALSE)
  sf <- with(fit$scale, c(sd_X / sd_Y, sd_X * sd_Z / sd_Y, sd_Z / sd_Y))
  expect_equal(ci$median, ci_raw$median * sf)
  ci90 <- credible_intervals(fit, level = 0.9)
  expect_true(all(ci90$lower >= ci$lower))
  expect_true(all(ci90$upper <= ci$upper))
  expect_error(credible_intervals(fit, level = 1), "'level' must lie")
})

test_that("sign_reversal uses paired draws on the standardised scale", {
  sr <- sign_reversal(fit, Z_observed = sim$Z)
  expect_true(sr$n_draws > 0)
  expect_length(sr$tau_draws, sr$n_draws)
  expect_true(sr$tau_lower <= sr$tau_median)
  expect_true(sr$tau_median <= sr$tau_upper)
  expect_true(sr$p_in_range >= 0 && sr$p_in_range <= 1)
  # Recompute from the raw draws: tau_std = (-beta_X / beta_XZ) / sd_Z
  ok <- abs(fit$draws$beta_XZ) >= 1e-8
  tau_manual <- (-fit$draws$beta_X[ok] / fit$draws$beta_XZ[ok]) /
    fit$scale$sd_Z
  expect_equal(sr$tau_draws, tau_manual)
  # Without Z_observed the stored range is used.
  sr2 <- sign_reversal(fit)
  expect_equal(sr2$Z_range_std, sr$Z_range_std)
})

test_that("the Gelman-Rubin helper behaves as specified", {
  set.seed(7)
  x <- rnorm(500)
  # With B = 0 the statistic equals sqrt((n - 1) / n), which is 1 up to
  # the finite-sample correction.
  rh <- gelman_rubin(list(x, x, x, x))
  expect_equal(rh, sqrt(499 / 500))
  expect_equal(rh, 1, tolerance = 5e-3)
  shifted <- list(rnorm(500, 0), rnorm(500, 5), rnorm(500, 10), rnorm(500, 15))
  expect_true(gelman_rubin(shifted) > 1)
  expect_true(is.na(gelman_rubin(list(x))))
  # Constant chains at different levels: total non-convergence.
  expect_equal(gelman_rubin(list(rep(0, 50), rep(1, 50))), Inf)
  # Constant and identical: undefined.
  expect_true(is.na(gelman_rubin(list(rep(1, 50), rep(1, 50)))))
})

test_that("batch-means MCSE is pooled in quadrature across chains", {
  set.seed(11)
  ch <- list(rbinom(600, 1, 0.5), rbinom(600, 1, 0.5))
  per <- lapply(ch, mcse_batch)
  pooled <- mcse_pooled(per)
  expect_equal(unname(pooled[["mcse"]]),
               sqrt(per[[1]][["mcse"]]^2 + per[[2]][["mcse"]]^2) / 2)
  expect_equal(unname(pooled[["ess"]]),
               per[[1]][["ess"]] + per[[2]][["ess"]])
  # Too short a chain gives NA rather than an error.
  expect_true(is.na(mcse_batch(rbinom(10, 1, 0.5))[["mcse"]]))
})

test_that("simulate_mrccc returns consistently shaped data", {
  s <- simulate_mrccc(n = 40, pG = 2, pH = 3, pV = 1, seed = 5)
  expect_equal(dim(s$X), c(40L, 1L))
  expect_equal(dim(s$G), c(40L, 2L))
  expect_equal(dim(s$H), c(40L, 3L))
  expect_equal(dim(s$V), c(40L, 1L))
  expect_equal(s$truth$tau, -0.3 / 0.3)
  s0 <- simulate_mrccc(n = 40, beta_XZ = 0, seed = 5)
  expect_true(is.na(s0$truth$tau))
  expect_error(simulate_mrccc(n = 0), "'n' must be a positive integer")
})
