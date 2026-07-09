# Confirms the reproducibility claim: a fixed seed produces identical NMF
# results run-to-run and across sequential vs parallel execution backends.
# These run locally / in a full environment; they are skipped without NMF.

small_matrix <- function() {
  data("example_expression_data", package = "NMFprofileR")
  m <- as.matrix(example_expression_data)
  v <- apply(m, 1, stats::var)
  m[order(v, decreasing = TRUE)[seq_len(min(30, nrow(m)))], , drop = FALSE]
}

test_that("NMF gives identical results run-to-run for a fixed seed", {
  skip_on_cran()
  skip_if_not_installed("NMF")

  m <- small_matrix()
  a <- NMF::nmf(m, rank = 2, method = "brunet", nrun = 2, seed = 123456,
                .options = list(parallel = 0))
  b <- NMF::nmf(m, rank = 2, method = "brunet", nrun = 2, seed = 123456,
                .options = list(parallel = 0))

  expect_equal(NMF::basis(a), NMF::basis(b))
  expect_equal(NMF::coef(a), NMF::coef(b))
  expect_equal(NMF::consensus(a), NMF::consensus(b))
})

test_that("the parallel backend does not change the seeded result", {
  skip_on_cran()
  skip_if_not_installed("NMF")
  if (parallel::detectCores() < 2L) skip("needs >= 2 cores")

  m <- small_matrix()
  # Use the package's own parallel path (nmf_fit), which loads NMF on the
  # workers so it works from within the package namespace.
  seq_fit <- nmf_fit(m, rank = 2, nrun = 2, seed = 123456, nmf_parallel = FALSE)
  par_fit <- nmf_fit(m, rank = 2, nrun = 2, seed = 123456, nmf_parallel = TRUE, n_cores = 2)

  expect_false(is.null(par_fit))
  expect_equal(NMF::basis(seq_fit), NMF::basis(par_fit))
  expect_equal(NMF::coef(seq_fit), NMF::coef(par_fit))
})