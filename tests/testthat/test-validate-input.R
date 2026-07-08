# Unit tests for the input-validation and duplicate-collapse helpers.

test_that("collapse_duplicate_genes collapses by max and preserves order", {
  m <- matrix(
    c(1, 2,
      9, 1,
      5, 5),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("A", "A", "B"), c("S1", "S2"))
  )
  out <- NMFprofileR:::collapse_duplicate_genes(m, method = "max")

  expect_equal(rownames(out), c("A", "B"))          # first-appearance order
  expect_equal(unname(out["A", ]), c(9, 2))         # per-sample maximum
  expect_equal(unname(out["B", ]), c(5, 5))
})

test_that("collapse_duplicate_genes can average duplicates", {
  m <- matrix(c(2, 4, 6, 8), nrow = 2, byrow = TRUE,
              dimnames = list(c("A", "A"), c("S1", "S2")))
  out <- NMFprofileR:::collapse_duplicate_genes(m, method = "mean")
  expect_equal(unname(out["A", ]), c(4, 6))
})

test_that("validate_expression_input requires gene-symbol rownames", {
  m <- matrix(1:6, nrow = 2)  # no rownames
  expect_error(NMFprofileR:::validate_expression_input(m), "rownames")
})

test_that("validate_expression_input errors on non-numeric input", {
  df <- data.frame(S1 = c("a", "b"), S2 = c("c", "d"), row.names = c("g1", "g2"))
  expect_error(NMFprofileR:::validate_expression_input(df), "numeric")
})

test_that("validate_expression_input collapses duplicates by max by default", {
  m <- matrix(c(1, 1, 9, 9, 3, 3), nrow = 3, byrow = TRUE,
              dimnames = list(c("A", "A", "B"), c("S1", "S2")))
  out <- NMFprofileR:::validate_expression_input(m)  # default collapse_max
  expect_equal(nrow(out), 2)
  expect_equal(unname(out["A", ]), c(9, 9))
})

test_that("validate_expression_input can error on duplicates instead", {
  m <- matrix(c(1, 1, 9, 9), nrow = 2, byrow = TRUE,
              dimnames = list(c("A", "A"), c("S1", "S2")))
  expect_error(
    NMFprofileR:::validate_expression_input(m, on_duplicate_genes = "error"),
    "duplicated gene symbol"
  )
})

test_that("validate_expression_input drops all-zero gene rows", {
  m <- matrix(c(0, 0, 5, 7), nrow = 2, byrow = TRUE,
              dimnames = list(c("zero", "keep"), c("S1", "S2")))
  out <- NMFprofileR:::validate_expression_input(m)
  expect_equal(rownames(out), "keep")
})