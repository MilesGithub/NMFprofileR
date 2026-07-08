# End-to-end smoke test on the bundled example data.
#
# This exercises the full pipeline, including live g:Profiler network calls, so
# it is deliberately kept out of the CRAN-style CI check. It runs locally via
# devtools::test() (which sets NOT_CRAN=true) when NMF is installed and the
# machine is online.

test_that("NMFprofileR() runs end-to-end on the example data", {
  skip_on_cran()
  skip_if_not_installed("NMF")
  skip_if_offline()

  data("example_expression_data", package = "NMFprofileR")

  out_prefix <- file.path(tempdir(), "smoke_NMF")

  res <- NMFprofileR(
    expression_data      = example_expression_data,
    nmf_rank             = 2:3,
    output_prefix        = out_prefix,
    nmf_nrun             = 2,
    expression_threshold = 0,   # keep all genes regardless of scale
    variance_quantile    = 0,
    verbose              = FALSE
  )

  expect_type(res, "list")
  expect_named(
    res,
    c("consolidated_summary_df", "rank_metrics", "nmf_rds", "output_dirs", "runtime")
  )
  expect_s3_class(res$consolidated_summary_df, "data.frame")
  expect_gt(nrow(res$consolidated_summary_df), 0)
  expect_true(dir.exists(paste0(out_prefix, "_Results")))
})
