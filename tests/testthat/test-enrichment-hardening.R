# Unit tests for the g:Profiler hardening helpers (P7): the query hash key, the
# retry wrapper, and run_gost_query()'s oversized-query guard and cache hit.
# None of these touch the network.

test_that("nmf_hash_key is deterministic and order-invariant for gene vectors", {
  k1 <- NMFprofileR:::nmf_hash_key(c("A", "B", "C"), "bg", "GO:BP")
  k2 <- NMFprofileR:::nmf_hash_key(c("C", "A", "B"), "bg", "GO:BP")   # reordered
  k3 <- NMFprofileR:::nmf_hash_key(c("A", "B", "D"), "bg", "GO:BP")   # different genes

  expect_type(k1, "character")
  expect_length(k1, 1L)
  expect_identical(k1, k2)
  expect_false(identical(k1, k3))
})

test_that("with_retry retries on error then succeeds, and gives up after exhausting", {
  attempts <- 0L
  flaky <- function() {
    attempts <<- attempts + 1L
    if (attempts < 3L) stop("transient")
    "ok"
  }
  expect_equal(NMFprofileR:::with_retry(flaky, retries = 5, delay = 0), "ok")
  expect_equal(attempts, 3L)

  expect_error(
    NMFprofileR:::with_retry(function() stop("always"), retries = 2, delay = 0),
    "always"
  )
})

test_that("run_gost_query skips oversized queries without calling the service", {
  expect_warning(
    res <- NMFprofileR:::run_gost_query(
      query = letters, organism = "hsapiens", sources = "GO:BP",
      correction = "g_SCS", cutoff = 0.05, background_genes = letters,
      query_size = 100L, max_query_size = 10L,
      retries = 0, retry_delay = 0, label = "big"
    ),
    "max_query_size"
  )
  expect_null(res)
})

test_that("run_gost_query returns a cached result without hitting the network", {
  cache <- file.path(tempdir(), paste0("gcache_", as.integer(runif(1, 1, 1e6))))
  dir.create(cache, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(cache, recursive = TRUE), add = TRUE)

  query <- c("A", "B", "C", "D", "E")
  bg <- c("A", "B", "C", "D", "E", "F", "G")
  sources <- "GO:BP"; correction <- "g_SCS"; cutoff <- 0.05; organism <- "hsapiens"

  # Pre-populate the cache with the exact key run_gost_query() computes.
  key <- NMFprofileR:::nmf_hash_key(query, bg, sources, correction, cutoff, organism)
  fake <- list(result = data.frame(term_id = "GO:1", p_value = 0.01), meta = list(version = "cached"))
  saveRDS(fake, file.path(cache, paste0("gost_", key, ".rds")))

  res <- NMFprofileR:::run_gost_query(
    query = query, organism = organism, sources = sources, correction = correction,
    cutoff = cutoff, background_genes = bg, query_size = length(query),
    max_query_size = 10000, cache_dir = cache, retries = 0, retry_delay = 0,
    label = "cache-hit"
  )
  expect_equal(res$meta$version, "cached")
  expect_s3_class(res$result, "data.frame")
})