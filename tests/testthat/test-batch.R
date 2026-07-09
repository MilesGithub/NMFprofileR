# Tests for the batch driver (P6). NMFprofileR() is mocked so these run without
# NMF or a network; they exercise consolidation, Run_ID tagging, skip_existing,
# on_error, and the reserved-argument guard.

fake_profile <- function(run_id, n_factors = 2L) {
  df <- tibble::tibble(
    Run_ID = run_id,
    Rank = 2L,
    Factor = paste0("Factor_", seq_len(n_factors)),
    Num_Samples = seq_len(n_factors)
  )
  structure(list(consolidated_summary_df = df), class = "nmf_profile")
}

test_that("run_nmf_batch consolidates cohorts and tags each with its Run_ID", {
  out <- file.path(tempdir(), paste0("batch_", as.integer(runif(1, 1, 1e6))))
  on.exit(unlink(out, recursive = TRUE), add = TRUE)

  cohorts <- list(A = matrix(1, 2, 2), B = matrix(1, 2, 2))

  res <- testthat::with_mocked_bindings(
    run_nmf_batch(cohorts, output_dir = out, nmf_rank = 2, nmf_nrun = 1),
    NMFprofileR = function(expression_data, output_prefix, run_id, ...) fake_profile(run_id),
    .package = "NMFprofileR"
  )

  expect_equal(sort(unique(res$consolidated$Run_ID)), c("A", "B"))
  expect_equal(nrow(res$consolidated), 4L)         # 2 factors x 2 cohorts
  expect_equal(nrow(res$failures), 0L)
  expect_length(res$skipped, 0L)
  # The consolidated summary is written to disk.
  expect_true(file.exists(file.path(out, "Batch_Consolidated_Summary.tsv")))
})

test_that("run_nmf_batch continues past a failing cohort and records it", {
  out <- file.path(tempdir(), paste0("batch_", as.integer(runif(1, 1, 1e6))))
  on.exit(unlink(out, recursive = TRUE), add = TRUE)

  cohorts <- list(good = matrix(1, 2, 2), bad = matrix(1, 2, 2))

  mock <- function(expression_data, output_prefix, run_id, ...) {
    if (identical(run_id, "bad")) stop("synthetic cohort failure")
    fake_profile(run_id)
  }

  res <- testthat::with_mocked_bindings(
    run_nmf_batch(cohorts, output_dir = out, on_error = "continue"),
    NMFprofileR = mock,
    .package = "NMFprofileR"
  )

  expect_equal(res$failures$Cohort, "bad")
  expect_match(res$failures$Reason, "synthetic cohort failure")
  expect_equal(unique(res$consolidated$Run_ID), "good")

  # on_error = "stop" propagates the first cohort error.
  expect_error(
    testthat::with_mocked_bindings(
      run_nmf_batch(cohorts, output_dir = out, on_error = "stop", skip_existing = FALSE),
      NMFprofileR = mock,
      .package = "NMFprofileR"
    ),
    "failed"
  )
})

test_that("run_nmf_batch skips cohorts whose results already exist", {
  out <- file.path(tempdir(), paste0("batch_", as.integer(runif(1, 1, 1e6))))
  on.exit(unlink(out, recursive = TRUE), add = TRUE)

  # Pre-create cohort A's results on disk (no Run_ID column, as an older run).
  a_dir <- file.path(out, "A_Results", "Summaries")
  dir.create(a_dir, showWarnings = FALSE, recursive = TRUE)
  readr::write_tsv(
    tibble::tibble(Rank = 2L, Factor = c("Factor_1", "Factor_2"), Num_Samples = c(3L, 1L)),
    file.path(a_dir, "Consolidated_Summary.tsv")
  )

  cohorts <- list(A = matrix(1, 2, 2), B = matrix(1, 2, 2))
  called <- character(0)
  mock <- function(expression_data, output_prefix, run_id, ...) {
    called <<- c(called, run_id)
    fake_profile(run_id)
  }

  res <- testthat::with_mocked_bindings(
    run_nmf_batch(cohorts, output_dir = out, skip_existing = TRUE),
    NMFprofileR = mock,
    .package = "NMFprofileR"
  )

  expect_equal(called, "B")                 # only B was actually run
  expect_equal(res$skipped, "A")
  # A is loaded from disk and tagged with its Run_ID; both cohorts present.
  expect_setequal(unique(res$consolidated$Run_ID), c("A", "B"))
})

test_that("run_nmf_batch rejects driver-managed arguments and bad cohort names", {
  expect_error(
    run_nmf_batch(list(A = matrix(1, 2, 2)), output_dir = tempdir(),
                  expression_data = matrix(1, 2, 2)),
    "do not pass|Do not pass"
  )
  expect_error(
    run_nmf_batch(list(matrix(1, 2, 2)), output_dir = tempdir()),
    "unique, non-empty names"
  )
})