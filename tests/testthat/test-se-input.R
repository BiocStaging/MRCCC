sim <- simulate_mrccc(n = 40, seed = 7)

m <- rbind(LIG1 = as.numeric(sim$X), REC1 = as.numeric(sim$Z))
colnames(m) <- paste0("donor", seq_len(ncol(m)))
se <- SummarizedExperiment::SummarizedExperiment(assays = list(logcounts = m))

test_that("se_row returns a named numeric donor vector", {
  x <- se_row(se, "LIG1")
  expect_true(is.numeric(x))
  expect_length(x, 40L)
  expect_identical(names(x), colnames(m))
  expect_equal(unname(x), as.numeric(sim$X))
})

test_that("se_row selects an assay by name", {
  expect_identical(se_row(se, "REC1", assay = "logcounts"),
                   se_row(se, "REC1"))
})

test_that("se_row rejects invalid input", {
  expect_error(se_row(se, "MISSING"), "not found")
  expect_error(se_row(m, "LIG1"), "must be a SummarizedExperiment")
  expect_error(se_row(se, c("LIG1", "REC1")), "single character")
  expect_error(se_row(se, "LIG1", assay = "counts"), "name of an assay")
  expect_error(se_row(se, "LIG1", assay = 3), "between 1 and 1")
})

test_that("a single-row SummarizedExperiment is accepted by mr_ccc", {
  se1 <- se[1L, ]
  fit <- suppressWarnings(
    mr_ccc(se1, sim$Z, sim$Y, sim$G, sim$H, sim$V,
           n_iter = 400, burn_in = 100, seed = 1)
  )
  expect_s3_class(fit, "mrccc_fit")
})

test_that("a multi-row SummarizedExperiment as X points to se_row", {
  expect_error(
    mr_ccc(se, sim$Z, sim$Y, sim$G, sim$H, sim$V,
           n_iter = 400, burn_in = 100),
    "se_row"
  )
})
