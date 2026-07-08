# write_files = FALSE must compute in memory and touch no disk.

test_that("setup_directories(create = FALSE) returns paths but creates nothing", {
  sandbox <- file.path(tempdir(), paste0("nodir_", as.integer(runif(1, 1, 1e6))))
  on.exit(unlink(sandbox, recursive = TRUE), add = TRUE)
  prefix <- file.path(sandbox, "Run")

  dirs <- NMFprofileR:::setup_directories(prefix, create = FALSE)

  expect_type(dirs, "list")
  expect_true(all(c("main", "plots", "summaries") %in% names(dirs)))
  expect_false(dir.exists(dirs$main))   # nothing created
  expect_false(dir.exists(sandbox))
})

test_that("NMFprofileR(write_files = FALSE) writes nothing and returns results", {
  skip_on_cran()
  skip_if_not_installed("NMF")
  skip_if_offline()

  data("example_expression_data", package = "NMFprofileR")

  sandbox <- file.path(tempdir(), paste0("nofiles_", as.integer(runif(1, 1, 1e6))))
  dir.create(sandbox, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(sandbox, recursive = TRUE), add = TRUE)

  res <- NMFprofileR(
    expression_data      = example_expression_data,
    nmf_rank             = 2,
    output_prefix        = file.path(sandbox, "NMF"),
    nmf_nrun             = 2,
    expression_threshold = 0,
    variance_quantile    = 0,
    write_files          = FALSE
  )

  expect_type(res, "list")
  expect_null(res$output_dirs)
  expect_s3_class(res$consolidated_summary_df, "data.frame")
  # The sandbox must remain empty: no *_Results tree, no log, no files at all.
  written <- list.files(sandbox, recursive = TRUE, all.files = TRUE, no.. = TRUE)
  expect_equal(length(written), 0L)
})