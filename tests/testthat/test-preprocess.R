# Unit tests for the internal preprocess_matrix() helper.
# Reached via the ::: operator since the helper is not exported.

test_that("non-finite values are imputed with 0", {
  m <- matrix(
    c(1, 2, NA, 4,
      5, 6, Inf, 8),
    nrow = 2, byrow = TRUE,
    dimnames = list(c("geneA", "geneB"), paste0("S", 1:4))
  )

  out <- NMFprofileR:::preprocess_matrix(m, threshold = 0, v_quantile = 0, filter_file = NULL)

  expect_false(any(!is.finite(out)))
  # The imputed positions should now be 0.
  expect_equal(out["geneA", "S3"], 0)
})

test_that("genes below the mean-expression threshold are removed", {
  m <- matrix(
    c(1, 1, 1, 1,          # geneLow  mean 1
      100, 100, 100, 100), # geneHigh mean 100
    nrow = 2, byrow = TRUE,
    dimnames = list(c("geneLow", "geneHigh"), paste0("S", 1:4))
  )

  out <- NMFprofileR:::preprocess_matrix(m, threshold = 10, v_quantile = 0, filter_file = NULL)

  expect_equal(rownames(out), "geneHigh")
})

test_that("low-variance genes are removed by the variance quantile filter", {
  m <- matrix(
    c(10, 10, 10, 10,   # geneFlat  variance 0
      1,  2,  3,  4,    # geneMid   modest variance
      0, 20,  0, 20),   # geneWide  high variance
    nrow = 3, byrow = TRUE,
    dimnames = list(c("geneFlat", "geneMid", "geneWide"), paste0("S", 1:4))
  )

  out <- NMFprofileR:::preprocess_matrix(m, threshold = 0, v_quantile = 0.5, filter_file = NULL)

  expect_false("geneFlat" %in% rownames(out)) # zero-variance gene dropped
  expect_true("geneWide" %in% rownames(out))  # highest-variance gene retained
})

test_that("negative values are clamped to 0", {
  m <- matrix(
    c(-1, 5, 5, 5,
      6, -2, 7, 8),
    nrow = 2, byrow = TRUE,
    dimnames = list(c("geneA", "geneB"), paste0("S", 1:4))
  )

  out <- NMFprofileR:::preprocess_matrix(m, threshold = 0, v_quantile = 0, filter_file = NULL)

  expect_true(all(out >= 0))
})

test_that("columns are preserved and output is non-negative", {
  m <- matrix(
    c(1, 2, 3, 4,
      5, 6, 7, 8,
      2, 4, 6, 8),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("geneA", "geneB", "geneC"), paste0("S", 1:4))
  )

  out <- NMFprofileR:::preprocess_matrix(m, threshold = 0, v_quantile = 0, filter_file = NULL)

  expect_equal(ncol(out), 4)
  expect_equal(colnames(out), paste0("S", 1:4))
  expect_true(all(out >= 0))
})

test_that("an external gene-filter file restricts the retained genes", {
  m <- matrix(
    c(1, 9, 1, 9,     # geneA
      2, 8, 2, 8,     # geneB
      3, 7, 3, 7),    # geneC
    nrow = 3, byrow = TRUE,
    dimnames = list(c("geneA", "geneB", "geneC"), paste0("S", 1:4))
  )

  filter_path <- tempfile(fileext = ".txt")
  writeLines(c("geneA", "geneC"), filter_path)
  on.exit(unlink(filter_path), add = TRUE)

  out <- NMFprofileR:::preprocess_matrix(m, threshold = 0, v_quantile = 0, filter_file = filter_path)

  expect_true(all(rownames(out) %in% c("geneA", "geneC")))
  expect_false("geneB" %in% rownames(out))
})
