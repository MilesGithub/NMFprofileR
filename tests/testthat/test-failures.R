# Failed ranks (P5) must be recorded in $failures and reported, not silently
# dropped. We force the per-rank fit to fail by mocking run_nmf_for_rank(), so
# the test is deterministic and needs no network (no rank succeeds, so no
# enrichment call is made).

test_that("failed ranks are captured in $failures and surfaced by print()", {
  skip_if_not_installed("NMF")

  set.seed(1)
  m <- matrix(abs(rnorm(8 * 12)) + 1, nrow = 8, ncol = 12,
              dimnames = list(paste0("G", 1:8), paste0("S", 1:12)))

  sandbox <- file.path(tempdir(), paste0("fail_", as.integer(runif(1, 1, 1e6))))
  on.exit(unlink(sandbox, recursive = TRUE), add = TRUE)

  fake_fit <- function(k, ...) {
    warning(sprintf("NMF fit failed for rank %s: synthetic failure", k), call. = FALSE)
    NULL
  }

  res <- testthat::with_mocked_bindings(
    suppressWarnings(
      NMFprofileR(
        expression_data      = m,
        nmf_rank             = 2:3,
        output_prefix        = file.path(sandbox, "Run"),
        nmf_nrun             = 1,
        expression_threshold = 0,
        variance_quantile    = 0,
        write_files          = FALSE
      )
    ),
    run_nmf_for_rank = fake_fit,
    .package = "NMFprofileR"
  )

  expect_s3_class(res$failures, "data.frame")
  expect_equal(res$failures$Rank, c(2L, 3L))
  expect_true(all(grepl("synthetic failure", res$failures$Reason)))
  # No rank succeeded, so there is no consolidated summary.
  expect_equal(nrow(res$consolidated_summary_df), 0L)

  out <- paste(utils::capture.output(print(res)), collapse = "\n")
  expect_match(out, "Ranks failed : 2, 3")
})