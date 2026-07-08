# Unit tests for the reproducibility-provenance helpers.

test_that("extract_gprofiler_version pulls the version from gost meta", {
  g1 <- list(result = data.frame(), meta = list(version = "e111_eg58_p18"))

  expect_equal(NMFprofileR:::extract_gprofiler_version(list(g1)), "e111_eg58_p18")
  # NULL entries are skipped
  expect_equal(NMFprofileR:::extract_gprofiler_version(list(NULL, g1)), "e111_eg58_p18")
  # a bare gost object (not wrapped in a list) is handled
  expect_equal(NMFprofileR:::extract_gprofiler_version(g1), "e111_eg58_p18")
  # nothing usable -> NA
  expect_true(is.na(NMFprofileR:::extract_gprofiler_version(list(NULL))))
  expect_true(is.na(NMFprofileR:::extract_gprofiler_version(NULL)))
})

test_that("write_run_provenance writes a manifest and session info", {
  tmp <- file.path(tempdir(), paste0("prov_", as.integer(runif(1, 1, 1e6))))
  dirs <- list(summaries = file.path(tmp, "Summaries"))
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  params <- list(
    nmf_seed          = 123456,
    nmf_method        = "brunet",
    gprofiler_sources = c("GO:BP", "REAC")
  )

  out <- NMFprofileR:::write_run_provenance(
    dirs, "Test", params, "e111_eg58_p18", as.difftime(3, units = "mins")
  )

  expect_true(file.exists(out$manifest))
  expect_true(file.exists(out$session_info))

  manifest <- readLines(out$manifest)
  expect_true(any(grepl("e111_eg58_p18", manifest)))        # captured version
  expect_true(any(grepl("123456", manifest)))               # a parameter value
  expect_true(any(grepl("GO:BP, REAC", manifest, fixed = TRUE))) # collapsed vector
})

test_that("write_run_provenance records 'not captured' when the version is NA", {
  tmp <- file.path(tempdir(), paste0("prov_", as.integer(runif(1, 1, 1e6))))
  dirs <- list(summaries = file.path(tmp, "Summaries"))
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  out <- NMFprofileR:::write_run_provenance(
    dirs, "Test", list(nmf_seed = 1), NA_character_, as.difftime(1, units = "secs")
  )
  manifest <- readLines(out$manifest)
  expect_true(any(grepl("not captured", manifest)))
})