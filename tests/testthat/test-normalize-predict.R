# Unit tests for the internal normalize_predict() helper, which coerces the
# varied return shapes of NMF::predict() into an integer vector.

test_that("an atomic vector of the expected length is returned as integers", {
  out <- NMFprofileR:::normalize_predict(c(1, 2, 1), expected_length = 3)
  expect_identical(out, c(1L, 2L, 1L))
})

test_that("a list with a $predict element is used directly", {
  out <- NMFprofileR:::normalize_predict(list(predict = c(2, 1, 3)), expected_length = 3)
  expect_identical(out, c(2L, 1L, 3L))
})

test_that("a features x k membership matrix is argmaxed row-wise", {
  m <- matrix(
    c(0.1, 0.9,
      0.8, 0.2,
      0.3, 0.7),
    nrow = 3, byrow = TRUE
  )
  out <- NMFprofileR:::normalize_predict(m, expected_length = 3)
  expect_identical(out, c(2L, 1L, 2L))
})

test_that("a transposed k x features membership matrix is argmaxed column-wise", {
  m <- matrix(
    c(0.1, 0.8, 0.3,
      0.9, 0.2, 0.7),
    nrow = 2, byrow = TRUE
  )
  out <- NMFprofileR:::normalize_predict(m, expected_length = 3)
  expect_identical(out, c(2L, 1L, 2L))
})

test_that("unrecognised input degrades to all-NA of the expected length", {
  out <- NMFprofileR:::normalize_predict(matrix(1:4, nrow = 2), expected_length = 3)
  expect_length(out, 3)
  expect_true(all(is.na(out)))
})

test_that("predictions shorter than expected are padded with NA", {
  out <- NMFprofileR:::normalize_predict(list(predict = c(1, 2)), expected_length = 4)
  expect_identical(out, c(1L, 2L, NA_integer_, NA_integer_))
})
