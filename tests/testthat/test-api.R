# Tests for the exported composable pipeline functions.

test_that("nmf_preprocess validates and filters end to end", {
  data("example_expression_data", package = "NMFprofileR")
  m <- nmf_preprocess(example_expression_data, expression_threshold = 0, variance_quantile = 0)
  expect_true(is.matrix(m))
  expect_gt(nrow(m), 2)
  expect_equal(ncol(m), ncol(example_expression_data))
  expect_true(all(m >= 0))
})

test_that("nmf_preprocess errors when nothing survives filtering", {
  data("example_expression_data", package = "NMFprofileR")
  expect_error(
    nmf_preprocess(example_expression_data, expression_threshold = 1e9),
    "Fewer than 2 genes"
  )
})

test_that("the composable stages round-trip a fit", {
  skip_on_cran()
  skip_if_not_installed("NMF")

  data("example_expression_data", package = "NMFprofileR")
  m <- nmf_preprocess(example_expression_data, expression_threshold = 0, variance_quantile = 0)
  fit <- nmf_fit(m, rank = 2, nrun = 2)
  skip_if(is.null(fit), "NMF fit unavailable")

  bg <- nmf_basis_genes(fit)
  expect_true(all(c("Gene", "Factor") %in% names(bg)))
  expect_true(all(as.character(bg$Factor) %in% paste0("Factor_", 1:2)))

  sa <- nmf_sample_assignments(fit)
  expect_equal(nrow(sa), ncol(m))
  expect_true(all(c("SampleID", "Dominant_Factor", "Silhouette_NMF") %in% names(sa)))
})

test_that("nmf_profile print and summary methods work", {
  obj <- structure(
    list(
      consolidated_summary_df = data.frame(Rank = 2L, Factor = "Factor_1"),
      fits = list("2" = NULL),
      sample_assignments = list("2" = data.frame(SampleID = c("a", "b"))),
      output_dirs = NULL,
      runtime = as.difftime(1, units = "mins"),
      provenance = list(gprofiler_version = "e111_test")
    ),
    class = "nmf_profile"
  )

  out <- capture.output(print(obj))
  expect_true(any(grepl("nmf_profile", out)))
  expect_true(any(grepl("in memory only", out)))   # output_dirs is NULL
  expect_true(any(grepl("e111_test", out)))
  expect_identical(summary(obj), obj$consolidated_summary_df)
})