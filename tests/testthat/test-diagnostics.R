# Tests for the rank-diagnostics helper.

test_that("nmf_rank_diagnostics returns tidy long metrics", {
  rm <- data.frame(
    rank = 2:4,
    cophenetic = c(0.99, 0.95, 0.88),
    dispersion = c(0.90, 0.80, 0.70),
    silhouette = c(0.60, 0.50, 0.40)
  )
  d <- nmf_rank_diagnostics(rm)

  expect_true(all(c("rank", "metric", "value") %in% names(d)))
  expect_setequal(unique(d$metric), c("cophenetic", "dispersion", "silhouette"))
  expect_equal(nrow(d), 9L)  # 3 metrics x 3 ranks
  expect_equal(d$value[d$metric == "cophenetic" & d$rank == 2], 0.99)
})

test_that("nmf_rank_diagnostics matches by prefix and skips missing metrics", {
  rm <- data.frame(rank = 2:3, cophenetic = c(0.9, 0.8), silhouette.consensus = c(0.5, 0.4))
  d <- nmf_rank_diagnostics(rm, metrics = c("cophenetic", "dispersion", "silhouette"))
  expect_setequal(unique(d$metric), c("cophenetic", "silhouette"))  # dispersion absent
  expect_equal(nrow(d), 4L)
})

test_that("nmf_rank_diagnostics handles empty input", {
  expect_equal(nrow(nmf_rank_diagnostics(data.frame())), 0L)
})