# Unit test for compute_sample_silhouette()'s graceful degradation when the
# NMF fit is unusable (no recoverable consensus matrix).

test_that("an unusable fit yields an NA silhouette frame keyed by colnames(H)", {
  H <- matrix(
    1:6, nrow = 2,
    dimnames = list(NULL, c("S1", "S2", "S3"))
  )

  # A NULL fit cannot produce a consensus matrix; the helper should catch this
  # internally and return NA silhouettes rather than error.
  out <- NMFprofileR:::compute_sample_silhouette(NULL, H, k = 2, verbose = FALSE)

  expect_s3_class(out, "data.frame")
  expect_equal(out$SampleID, c("S1", "S2", "S3"))
  expect_true(all(is.na(out$Silhouette_NMF)))
})
