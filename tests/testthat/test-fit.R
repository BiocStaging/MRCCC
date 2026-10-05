sim <- simulate_mrccc(n = 150, seed = 1)

fit <- suppressWarnings(
  mr_ccc(sim$X, sim$Z, sim$Y, sim$G, sim$H, sim$V,
         n_iter = 1500, burn_in = 300, seed = 1)
)

test_that("mr_ccc returns a well-formed mrccc_fit object", {
  expect_s3_class(fit, "mrccc_fit")
  expect_named(fit, c("pip", "pip_mcse", "pip_ess", "estimates", "draws",
                      "scale", "diagnostics", "instruments", "settings",
                      "call"))
  expect_true(is.numeric(fit$pip) && length(fit$pip) == 1L)
  expect_true(fit$pip >= 0 && fit$pip <= 1)
  expect_true(is.na(fit$pip_mcse) || fit$pip_mcse >= 0)
})

test_that("the default is four chains", {
  expect_equal(fit$settings$n_chains, 4L)
})

test_that("draws are pooled with the expected length", {
  n_keep <- fit$settings$n_keep
  expect_equal(n_keep, 1200L)
  for (nm in c("beta_X", "beta_XZ", "beta_Z", "gamma", "mu", "loglik")) {
    expect_length(fit$draws[[nm]], n_keep * fit$settings$n_chains)
    expect_true(all(is.finite(fit$draws[[nm]])))
  }
  expect_true(all(fit$draws$gamma %in% c(0, 1)))
  expect_equal(mean(fit$draws$gamma), fit$pip)
})

test_that("estimates and diagnostics tables have the expected shape", {
  expect_s3_class(fit$estimates, "data.frame")
  expect_equal(nrow(fit$estimates), 3L)
  expect_equal(fit$estimates$parameter, c("beta_X", "beta_XZ", "beta_Z"))
  expect_named(fit$estimates, c("parameter", "raw_mean", "std_mean"))
  sf <- with(fit$scale, c(sd_X / sd_Y, sd_X * sd_Z / sd_Y, sd_Z / sd_Y))
  expect_equal(fit$estimates$std_mean, fit$estimates$raw_mean * sf)

  expect_s3_class(fit$diagnostics, "data.frame")
  expect_equal(fit$diagnostics$parameter,
               c("beta_X", "beta_XZ", "gamma", "loglik"))
})

test_that("first-stage F is reported and is large for strong simulated instruments", {
  fi <- fit$instruments
  expect_s3_class(fi, "data.frame")
  expect_equal(fi$exposure, c("X", "Z"))
  expect_equal(fi$n_instruments, c(ncol(sim$G), ncol(sim$H)))
  expect_true(all(is.finite(fi$F_statistic)))
  expect_true(all(fi$F_statistic > 10))
  expect_false(any(fi$weak))

  # Weak instruments: pure-noise genotypes give F near 1 and a message.
  set.seed(4)
  G_noise <- matrix(stats::rnorm(nrow(sim$G) * ncol(sim$G)), nrow(sim$G))
  expect_message(
    fit_w <- suppressWarnings(
      mr_ccc(sim$X, sim$Z, sim$Y, G_noise, sim$H, sim$V,
             n_iter = 400, burn_in = 100, n_chains = 1, seed = 1)),
    "First-stage F below 10 for X"
  )
  expect_true(fit_w$instruments$weak[fit_w$instruments$exposure == "X"])
  expect_false(fit_w$instruments$weak[fit_w$instruments$exposure == "Z"])
})

test_that("R-hat is NA for a single chain", {
  fit1 <- suppressWarnings(
    mr_ccc(sim$X, sim$Z, sim$Y, sim$G, sim$H, sim$V,
           n_iter = 1000, burn_in = 200, n_chains = 1, seed = 1)
  )
  expect_equal(fit1$settings$n_chains, 1L)
  expect_true(all(is.na(fit1$diagnostics$rhat)))
})

test_that("the log-likelihood is finite and not constant", {
  ll <- fit$draws$loglik
  expect_true(all(is.finite(ll)))
  expect_gt(stats::var(ll), 0)
})

test_that("a seed makes the fit reproducible", {
  fit2 <- suppressWarnings(
    mr_ccc(sim$X, sim$Z, sim$Y, sim$G, sim$H, sim$V,
           n_iter = 1500, burn_in = 300, seed = 1)
  )
  expect_identical(fit2$pip, fit$pip)
  expect_identical(fit2$draws$beta_XZ, fit$draws$beta_XZ)
})

test_that("shifting instrument or covariate columns leaves the fit unchanged", {
  G2 <- sim$G + 2
  V2 <- sim$V + 5
  fit_shift <- suppressWarnings(
    mr_ccc(sim$X, sim$Z, sim$Y, G2, sim$H, V2,
           n_iter = 1500, burn_in = 300, seed = 1)
  )
  expect_equal(fit_shift$pip, fit$pip, tolerance = 1e-8)
  expect_equal(fit_shift$draws$beta_XZ, fit$draws$beta_XZ, tolerance = 1e-6)
  expect_equal(fit_shift$settings$prior_scale, fit$settings$prior_scale,
               tolerance = 1e-8)
})

test_that("vector inputs are accepted and give the same result as matrices", {
  fit3 <- suppressWarnings(
    mr_ccc(as.numeric(sim$X), as.numeric(sim$Z), as.numeric(sim$Y),
           sim$G, sim$H, sim$V, n_iter = 1500, burn_in = 300, seed = 1)
  )
  expect_identical(fit3$pip, fit$pip)
})

test_that("multiple chains give finite R-hat and pooled draws", {
  fit4 <- suppressWarnings(
    mr_ccc(sim$X, sim$Z, sim$Y, sim$G, sim$H, sim$V,
           n_iter = 1000, burn_in = 200, n_chains = 2, seed = 2)
  )
  expect_equal(fit4$settings$n_chains, 2L)
  expect_length(fit4$draws$gamma, 2L * fit4$settings$n_keep)
  rh <- fit4$diagnostics$rhat
  expect_true(all(is.finite(rh[fit4$diagnostics$parameter != "gamma"])))
  expect_true(all(rh[fit4$diagnostics$parameter != "gamma"] > 0))
})

test_that("thinning reduces the number of retained draws", {
  fit5 <- suppressWarnings(
    mr_ccc(sim$X, sim$Z, sim$Y, sim$G, sim$H, sim$V,
           n_iter = 700, burn_in = 100, thin = 3, n_chains = 1, seed = 1)
  )
  expect_equal(fit5$settings$n_keep, 200L)
  expect_length(fit5$draws$beta_X, 200L)
})

test_that("print and summary methods run and summary returns a table", {
  expect_output(print(fit), "MR-CCC fit")
  expect_output(tab <- summary(fit), "Posterior inclusion probability")
  expect_s3_class(tab, "data.frame")
  expect_equal(tab$parameter, c("beta_X", "beta_XZ", "beta_Z", "pip"))
  expect_equal(tab$raw_mean[4], fit$pip)
})

test_that("plot methods return a ggplot object", {
  p1 <- plot_effect_curve(fit, Z_observed = sim$Z)
  p2 <- plot(fit)
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})

test_that("the compiled sampler honours n_keep and verbose", {
  ctr <- function(v) matrix(v - mean(v), ncol = 1)
  ctr_cols <- function(M) scale(M, center = TRUE, scale = FALSE)
  Gc <- ctr_cols(sim$G); Hc <- ctr_cols(sim$H); Vc <- ctr_cols(sim$V)
  out <- mr_ccc_gibbs(ctr(sim$X), ctr(sim$Z), ctr(sim$Y),
                      Gc, Hc, Vc,
                      n_iter = 300, burn_in = 100, thin = 2)
  expect_equal(out$n_keep, 100L)
  expect_length(out$gamma_draws, 100L)
  expect_length(out$loglik_draws, 100L)
  expect_true(all(is.finite(out$loglik_draws)))
  expect_silent(
    mr_ccc_gibbs(ctr(sim$X), ctr(sim$Z), ctr(sim$Y), Gc, Hc, Vc,
                 n_iter = 1000, burn_in = 100, thin = 10)
  )
  expect_output(
    mr_ccc_gibbs(ctr(sim$X), ctr(sim$Z), ctr(sim$Y), Gc, Hc, Vc,
                 n_iter = 1000, burn_in = 100, thin = 10, verbose = TRUE),
    "Iteration 1000 / 1000"
  )
  expect_named(out$prior_scale, c("d_X", "d_XZ", "d_Z"))
})

test_that("the compiled sampler rejects invalid arguments", {
  ctr <- function(v) matrix(v - mean(v), ncol = 1)
  ctr_cols <- function(M) scale(M, center = TRUE, scale = FALSE)
  Gc <- ctr_cols(sim$G); Hc <- ctr_cols(sim$H); Vc <- ctr_cols(sim$V)
  expect_error(
    mr_ccc_gibbs(ctr(sim$X), ctr(sim$Z), ctr(sim$Y), Gc, Hc, Vc,
                 n_iter = 300, burn_in = 100, thin = 0),
    "thin"
  )
  expect_error(
    mr_ccc_gibbs(ctr(sim$X), ctr(sim$Z), ctr(sim$Y), Gc, Hc, Vc,
                 n_iter = 100, burn_in = 100),
    "burn_in"
  )
  expect_error(
    mr_ccc_gibbs(ctr(sim$X), ctr(sim$Z), ctr(sim$Y), Gc[-1, , drop = FALSE],
                 Hc, Vc, n_iter = 300, burn_in = 100),
    "same number of rows"
  )
  expect_error(
    mr_ccc_gibbs(ctr(sim$X), ctr(sim$Z), ctr(sim$Y), Gc, Hc, Vc,
                 n_iter = 300, burn_in = 100, nu1 = 0),
    "nu1"
  )
})
