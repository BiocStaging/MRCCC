sim <- simulate_mrccc(n = 60, seed = 3)

# n_iter and burn_in are formal arguments here rather than fixed inside the
# body, so that a test can override either without the call matching the
# same argument twice.
small_fit <- function(..., n_iter = 50, burn_in = 10) {
  mr_ccc(sim$X, sim$Z, sim$Y, sim$G, sim$H, sim$V,
         n_iter = n_iter, burn_in = burn_in, ...)
}

test_that("mismatched row counts are rejected", {
  expect_error(
    mr_ccc(sim$X[-1, , drop = FALSE], sim$Z, sim$Y, sim$G, sim$H, sim$V,
           n_iter = 50, burn_in = 10),
    "same number of rows"
  )
  expect_error(
    mr_ccc(sim$X, sim$Z, sim$Y, sim$G[-1, ], sim$H, sim$V,
           n_iter = 50, burn_in = 10),
    "same number of rows"
  )
})

test_that("missing and non-finite values are rejected", {
  Xna <- sim$X; Xna[5, 1] <- NA
  expect_error(
    mr_ccc(Xna, sim$Z, sim$Y, sim$G, sim$H, sim$V, n_iter = 50, burn_in = 10),
    "'X' contains NA"
  )
  Ginf <- sim$G; Ginf[2, 1] <- Inf
  expect_error(
    mr_ccc(sim$X, sim$Z, sim$Y, Ginf, sim$H, sim$V, n_iter = 50, burn_in = 10),
    "'G' contains NA"
  )
})

test_that("non-numeric inputs are rejected", {
  expect_error(
    mr_ccc(as.character(sim$X), sim$Z, sim$Y, sim$G, sim$H, sim$V,
           n_iter = 50, burn_in = 10),
    "'X' must be numeric"
  )
  expect_error(
    mr_ccc(sim$X, sim$Z, sim$Y, matrix("a", 60, 2), sim$H, sim$V,
           n_iter = 50, burn_in = 10),
    "'G' must be a numeric matrix"
  )
})

test_that("multi-column X, Z or Y is rejected", {
  expect_error(
    mr_ccc(cbind(sim$X, sim$X), sim$Z, sim$Y, sim$G, sim$H, sim$V,
           n_iter = 50, burn_in = 10),
    "exactly one column"
  )
})

test_that("V = NULL is rejected with an informative message", {
  expect_error(
    mr_ccc(sim$X, sim$Z, sim$Y, sim$G, sim$H, V = NULL,
           n_iter = 50, burn_in = 10),
    "'V' must be a numeric matrix with at least one covariate column"
  )
})

test_that("too few donors are rejected", {
  idx <- seq_len(8)
  expect_error(
    mr_ccc(sim$X[idx, , drop = FALSE], sim$Z[idx, , drop = FALSE],
           sim$Y[idx, , drop = FALSE], sim$G[idx, ], sim$H[idx, ],
           sim$V[idx, ], n_iter = 50, burn_in = 10),
    "At least 10 donors"
  )
})

test_that("constant, duplicated and saturating design columns are rejected", {
  Gconst <- sim$G; Gconst[, 2] <- 1
  expect_error(
    mr_ccc(sim$X, sim$Z, sim$Y, Gconst, sim$H, sim$V,
           n_iter = 50, burn_in = 10),
    "'G' has constant column"
  )

  Vconst <- cbind(sim$V, 1)
  expect_error(
    mr_ccc(sim$X, sim$Z, sim$Y, sim$G, sim$H, Vconst,
           n_iter = 50, burn_in = 10),
    "already includes an intercept"
  )

  Gdup <- sim$G; Gdup[, 3] <- Gdup[, 1]
  expect_error(
    mr_ccc(sim$X, sim$Z, sim$Y, Gdup, sim$H, sim$V,
           n_iter = 50, burn_in = 10),
    "duplicated column"
  )

  # More instrument columns than donors saturates the first stage.
  set.seed(11)
  idx   <- seq_len(12)
  Gwide <- matrix(stats::rnorm(12 * 15), nrow = 12)
  expect_error(
    mr_ccc(sim$X[idx, , drop = FALSE], sim$Z[idx, , drop = FALSE],
           sim$Y[idx, , drop = FALSE], Gwide, sim$H[idx, ], sim$V[idx, ],
           n_iter = 50, burn_in = 10),
    "saturated"
  )
})

test_that("too few donors for the number of columns is rejected", {
  idx <- seq_len(16)
  expect_error(
    mr_ccc(sim$X[idx, , drop = FALSE], sim$Z[idx, , drop = FALSE],
           sim$Y[idx, , drop = FALSE], sim$G[idx, ], sim$H[idx, ],
           sim$V[idx, ], n_iter = 50, burn_in = 10),
    "Too few donors for the model"
  )
})

test_that("a configuration retaining fewer than two draws is rejected", {
  expect_error(small_fit(n_iter = 20, burn_in = 10, thin = 40),
               "retain only")
})

test_that("a non-integer seed is rejected", {
  expect_error(small_fit(seed = 1.5), "'seed' must be a whole number")
})

test_that("MCMC settings are validated", {
  expect_error(
    mr_ccc(sim$X, sim$Z, sim$Y, sim$G, sim$H, sim$V,
           n_iter = 100, burn_in = 100),
    "'n_iter' \\(100\\) must be greater than 'burn_in'"
  )
  expect_error(
    mr_ccc(sim$X, sim$Z, sim$Y, sim$G, sim$H, sim$V,
           n_iter = 100, burn_in = 200),
    "must be greater than 'burn_in'"
  )
  expect_error(small_fit(thin = 0), "'thin' must be an integer")
  expect_error(small_fit(thin = 1.5), "'thin' must be an integer")
  expect_error(small_fit(n_chains = 0), "'n_chains' must be an integer")
  expect_error(small_fit(init_scale = -1), "'init_scale' must be non-negative")
  expect_error(small_fit(burn_in = -1), "'burn_in' must be a non-negative")
})

test_that("pip_threshold must lie strictly inside (0, 1)", {
  expect_error(small_fit(pip_threshold = 0), "'pip_threshold' must lie")
  expect_error(small_fit(pip_threshold = 1), "'pip_threshold' must lie")
  expect_error(small_fit(pip_threshold = 1.5), "'pip_threshold' must lie")
})

test_that("hyperparameter constructor validates its inputs", {
  expect_s3_class(mrccc_hyperparameters(), "mrccc_hyper")
  expect_null(mrccc_hyperparameters()$g)
  expect_equal(mrccc_hyperparameters(g = 50)$g, 50)
  expect_error(mrccc_hyperparameters(nu1 = 1), "'nu1' must lie")
  expect_error(mrccc_hyperparameters(a_sigma = 0), "must be positive")
  expect_error(mrccc_hyperparameters(g = -1), "'g' must be positive")
  expect_error(mrccc_hyperparameters(ridge = -1e-8), "non-negative")
  expect_error(small_fit(hyper = list(a_sigma = 3)),
               "'hyper' must be created by mrccc_hyperparameters")
})

test_that("helper functions reject objects that are not fits", {
  expect_error(credible_intervals(list()), "must be an object of class")
  expect_error(sign_probabilities(1), "must be an object of class")
  expect_error(sign_reversal("a"), "must be an object of class")
  expect_error(plot_effect_curve(NULL), "must be an object of class")
})

test_that("rank-deficient designs are rejected after centring", {
  sim <- simulate_mrccc(n = 120, seed = 3)
  G_bad <- cbind(sim$G, 2 - sim$G[, 1])
  expect_error(
    mr_ccc(sim$X, sim$Z, sim$Y, G_bad, sim$H, sim$V,
           n_iter = 300, burn_in = 100, seed = 1),
    "rank deficient"
  )
})

test_that("a supplied seed leaves the caller's random number state intact", {
  sim <- simulate_mrccc(n = 120, seed = 3)
  set.seed(99)
  before <- runif(1)
  set.seed(99)
  invisible(suppressWarnings(suppressMessages(
    mr_ccc(sim$X, sim$Z, sim$Y, sim$G, sim$H, sim$V,
           n_iter = 300, burn_in = 100, seed = 1))))
  expect_identical(runif(1), before)
})
