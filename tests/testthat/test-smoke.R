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
    expression_threshold = 10,   # a small, fast gene set keeps the smoke test quick
    variance_quantile    = 0.99,
    run_id               = "smoke_run",
    verbose              = FALSE
  )

  expect_s3_class(res, "nmf_profile")
  expect_true(all(
    c("consolidated_summary_df", "rank_metrics", "fits", "basis_genes",
      "marker_genes", "sample_assignments", "enrichment", "failures",
      "provenance") %in% names(res)
  ))
  expect_true(all(c("per_factor", "combined", "markers") %in% names(res$enrichment)))
  expect_s3_class(res$consolidated_summary_df, "data.frame")
  expect_gt(nrow(res$consolidated_summary_df), 0)
  # in-memory objects are returned, keyed by rank
  expect_named(res$fits, c("2", "3"))
  expect_equal(nrow(res$sample_assignments[["2"]]), ncol(example_expression_data))
  expect_true(dir.exists(paste0(out_prefix, "_Results")))

  # Every rank fitted, so there are no failures (P5).
  expect_s3_class(res$failures, "data.frame")
  expect_equal(nrow(res$failures), 0L)

  # run_id is stamped as a leading column (P2).
  expect_true("Run_ID" %in% names(res$consolidated_summary_df))
  expect_equal(names(res$consolidated_summary_df)[1], "Run_ID")
  expect_true(all(res$consolidated_summary_df$Run_ID == "smoke_run"))

  # Markers are opt-in and off by default now (P3): keys exist but are empty.
  expect_equal(length(res$marker_genes), 0L)

  # Single-object bundle round-trips and the manifest lists it (P2).
  bundle <- file.path(paste0(out_prefix, "_Results"), "smoke_NMF_nmf_profile.rds")
  expect_true(file.exists(bundle))
  reloaded <- readRDS(bundle)
  expect_s3_class(reloaded, "nmf_profile")
  expect_equal(reloaded$consolidated_summary_df$Run_ID, res$consolidated_summary_df$Run_ID)

  manifest <- file.path(paste0(out_prefix, "_Results"), "Summaries", "manifest.tsv")
  expect_true(file.exists(manifest))
  man <- readr::read_tsv(manifest, show_col_types = FALSE)
  expect_true("profile_bundle" %in% man$type)
  expect_true("core_rds" %in% man$type)
})
