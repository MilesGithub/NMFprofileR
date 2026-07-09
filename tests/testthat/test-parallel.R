# The opt-in parallel NMF path (P8) must return exactly the same factorization
# as the sequential path for a fixed seed. Kept out of the CRAN-style check
# (spins up a cluster); runs locally and under coverage.

test_that("parallel nmf_fit reproduces the sequential fit for a fixed seed", {
  skip_on_cran()
  skip_if_not_installed("NMF")
  if (parallel::detectCores() < 2L) skip("needs >= 2 cores")

  set.seed(1)
  m <- matrix(abs(rnorm(40 * 16)) + 1, nrow = 40, ncol = 16,
              dimnames = list(paste0("G", 1:40), paste0("S", 1:16)))

  seq_fit <- nmf_fit(m, rank = 2, nrun = 3, seed = 42, nmf_parallel = FALSE)

  # Capture warnings so we can assert the parallel path actually ran rather than
  # falling back to sequential (which would make the equality check trivial).
  warns <- character(0)
  par_fit <- withCallingHandlers(
    nmf_fit(m, rank = 2, nrun = 3, seed = 42, nmf_parallel = TRUE, n_cores = 2),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_false(is.null(seq_fit))
  expect_false(is.null(par_fit))
  expect_false(any(grepl("falling back to sequential", warns)))
  expect_equal(NMF::basis(seq_fit), NMF::basis(par_fit))
  expect_equal(NMF::coef(seq_fit), NMF::coef(par_fit))
})

test_that("nmf_parallel with nrun = 1 falls back to a sequential fit", {
  skip_on_cran()
  skip_if_not_installed("NMF")

  set.seed(2)
  m <- matrix(abs(rnorm(30 * 12)) + 1, nrow = 30, ncol = 12,
              dimnames = list(paste0("G", 1:30), paste0("S", 1:12)))

  seq_fit <- nmf_fit(m, rank = 2, nrun = 1, seed = 7, nmf_parallel = FALSE)
  par_fit <- nmf_fit(m, rank = 2, nrun = 1, seed = 7, nmf_parallel = TRUE)

  expect_equal(NMF::basis(seq_fit), NMF::basis(par_fit))
})