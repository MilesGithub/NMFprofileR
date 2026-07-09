# Tests for the scientific tools (v0.5.0): align_factors(), nmf_project(),
# nmf_stability(), and the basis_gene_method switch. The alignment and
# projection tests are network-free and NMF-free (they use synthetic matrices).

make_basis <- function(seed = 1) {
  set.seed(seed)
  genes <- paste0("G", 1:50)
  matrix(abs(rnorm(50 * 3)), nrow = 50, ncol = 3,
         dimnames = list(genes, paste0("Factor_", 1:3)))
}

test_that("align_factors recovers a known factor permutation (greedy)", {
  W1 <- make_basis()
  perm <- c(3L, 1L, 2L)                       # W2 factor i == W1 factor perm[i]
  W2 <- W1[, perm] + matrix(rnorm(150, sd = 1e-3), 50, 3)
  colnames(W2) <- paste0("Factor_", 1:3)

  aln <- align_factors(list(A = W1, B = W2), method = "greedy", min_cor = 0.5)

  expect_equal(aln$reference, "A")
  expect_equal(nrow(aln$correlations), 9L)    # 3 x 3 pairwise
  expect_equal(nrow(aln$matches), 3L)

  m <- aln$matches
  expect_equal(m$ref_factor[m$factor == "Factor_1"], "Factor_3")
  expect_equal(m$ref_factor[m$factor == "Factor_2"], "Factor_1")
  expect_equal(m$ref_factor[m$factor == "Factor_3"], "Factor_2")
  expect_true(all(m$correlation > 0.99))
})

test_that("align_factors hungarian agrees with greedy on a clean permutation", {
  skip_if_not_installed("clue")
  W1 <- make_basis(2)
  perm <- c(2L, 3L, 1L)
  W2 <- W1[, perm]
  colnames(W2) <- paste0("Factor_", 1:3)

  h <- align_factors(list(A = W1, B = W2), method = "hungarian")$matches
  expect_equal(nrow(h), 3L)
  expect_equal(h$ref_factor[h$factor == "Factor_1"], "Factor_2")
  expect_equal(h$ref_factor[h$factor == "Factor_3"], "Factor_1")
})

test_that("align_factors validates its input", {
  W1 <- make_basis()
  expect_error(align_factors(list(W1)), "at least two")
  expect_error(align_factors(list(W1, W1)), "unique, non-empty names")
})

test_that("nmf_project scores new samples onto fixed factors", {
  W <- matrix(0, nrow = 20, ncol = 2,
              dimnames = list(paste0("G", 1:20), c("Factor_1", "Factor_2")))
  set.seed(1)
  W[1:10, 1] <- runif(10, 1, 2)
  W[11:20, 2] <- runif(10, 1, 2)

  V <- cbind(
    S1 = as.numeric(W %*% c(5, 0)),   # pure factor 1
    S2 = as.numeric(W %*% c(0, 5)),   # pure factor 2
    S3 = as.numeric(W %*% c(3, 1))    # mostly factor 1
  )
  rownames(V) <- rownames(W)

  proj <- nmf_project(W, V)

  expect_s3_class(proj, "data.frame")
  expect_equal(proj$SampleID, c("S1", "S2", "S3"))
  expect_true(all(c("Factor_1", "Factor_2", "Dominant_Factor") %in% names(proj)))
  expect_equal(proj$Dominant_Factor, c("Factor_1", "Factor_2", "Factor_1"))
  # coefficients are non-negative
  expect_true(all(proj$Factor_1 >= 0) && all(proj$Factor_2 >= 0))
})

test_that("nmf_project aligns genes and errors when too few are shared", {
  W <- matrix(runif(20), nrow = 10, ncol = 2,
              dimnames = list(paste0("G", 1:10), c("Factor_1", "Factor_2")))
  V <- matrix(runif(16), nrow = 8, ncol = 2,
              dimnames = list(paste0("G", 3:10), c("N1", "N2")))   # missing G1,G2
  # Missing basis genes are reported via a cli alert (not an R warning); the
  # projection still succeeds on the shared genes.
  proj <- nmf_project(W, V)
  expect_equal(proj$SampleID, c("N1", "N2"))
  expect_true("Dominant_Factor" %in% names(proj))

  V_bad <- matrix(runif(4), nrow = 2, ncol = 2,
                  dimnames = list(c("Z1", "G1"), c("N1", "N2")))    # 1 shared gene
  expect_error(nmf_project(W, V_bad), "Fewer than 2 shared")
})

test_that("nmf_stability returns a per-factor stability tibble", {
  skip_on_cran()
  skip_if_not_installed("NMF")

  set.seed(3)
  m <- matrix(abs(rnorm(40 * 20)) + 1, nrow = 40, ncol = 20,
              dimnames = list(paste0("G", 1:40), paste0("S", 1:20)))

  stab <- nmf_stability(m, rank = 2, nrun = 1, nboot = 3, subsample_frac = 0.8)

  expect_s3_class(stab, "tbl_df")
  expect_equal(stab$Factor, c("Factor_1", "Factor_2"))
  finite <- stab$Stability[is.finite(stab$Stability)]
  expect_true(all(finite >= 0 & finite <= 1.0000001))
})

test_that("basis_gene_method = 'specificity' uses marker genes as the primary output", {
  skip_on_cran()
  skip_if_not_installed("NMF")

  set.seed(4)
  m <- matrix(abs(rnorm(40 * 16)) + 1, nrow = 40, ncol = 16,
              dimnames = list(paste0("G", 1:40), paste0("S", 1:16)))
  sandbox <- file.path(tempdir(), paste0("bgm_", as.integer(runif(1, 1, 1e6))))
  on.exit(unlink(sandbox, recursive = TRUE), add = TRUE)

  # Mock enrichment so the pipeline needs no network.
  res <- testthat::with_mocked_bindings(
    NMFprofileR(
      expression_data = m, nmf_rank = 2, output_prefix = file.path(sandbox, "R"),
      nmf_nrun = 1, expression_threshold = 0, variance_quantile = 0,
      write_files = FALSE, basis_gene_method = "specificity"
    ),
    perform_enrichment = function(...) list(),
    perform_combined_enrichment = function(...) data.frame(),
    .package = "NMFprofileR"
  )

  # The primary basis genes must equal the Kim-Park markers of the same fit.
  expect_equal(
    sort(as.character(res$basis_genes[["2"]]$Gene)),
    sort(as.character(nmf_marker_genes(res$fits[["2"]])$Gene))
  )
})