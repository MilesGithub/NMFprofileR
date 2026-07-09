# Unit tests for the self-describing-output helpers (P2). These are pure
# file-bookkeeping functions, so they run without NMF or a network.

test_that("classify_output_file maps file names to the expected types", {
  cls <- NMFprofileR:::classify_output_file
  expect_equal(cls("Run_nmf_profile.rds"), "profile_bundle")
  expect_equal(cls("NMF_Result_Object_Rank_k3.rds"), "core_rds")
  expect_equal(cls("Marker_Genes_Rank_k3.tsv"), "marker_genes")
  expect_equal(cls("Basis_Genes_Rank_k3.tsv"), "basis_genes")
  expect_equal(cls("Sample_Assignments_Rank_k3.tsv"), "sample_assignments")
  expect_equal(cls("Enrichment_Rank_k3_Factor_1.tsv"), "enrichment")
  expect_equal(cls("Marker_Enrichment_Rank_k3_Factor_1.tsv"), "enrichment")
  expect_equal(cls("02_Factor_Summary_Plot_Rank_k3.pdf"), "plot")
  expect_equal(cls("Run_Session_Info.txt"), "provenance")
  expect_equal(cls("Run_Run_Manifest.txt"), "provenance")
  expect_equal(cls("Run_run_log.txt"), "log")
  expect_equal(cls("Consolidated_Summary.tsv"), "summary")
  expect_equal(cls("something_unexpected.dat"), "other")
})

test_that("write_output_manifest lists every produced file, relative and typed", {
  sandbox <- file.path(tempdir(), paste0("manifest_", as.integer(runif(1, 1, 1e6))))
  on.exit(unlink(sandbox, recursive = TRUE), add = TRUE)
  prefix <- file.path(sandbox, "Run")

  dirs <- NMFprofileR:::setup_directories(prefix, create = TRUE)

  # Drop representative dummy outputs into the tree.
  writeLines("x", file.path(dirs$nmf_core, "NMF_Result_Object_Rank_k2.rds"))
  writeLines("x", file.path(dirs$basis_genes, "Basis_Genes_Rank_k2.tsv"))
  writeLines("x", file.path(dirs$basis_genes, "Marker_Genes_Rank_k2.tsv"))
  writeLines("x", file.path(dirs$enrichment, "Enrichment_Rank_k2_Factor_1.tsv"))
  writeLines("x", file.path(dirs$sample_assignments, "Sample_Assignments_Rank_k2.tsv"))
  writeLines("x", file.path(dirs$plots, "02_Rank_Survey_Plot.pdf"))
  writeLines("x", file.path(dirs$summaries, "Consolidated_Summary.tsv"))
  writeLines("x", file.path(dirs$main, "Run_run_log.txt"))
  bundle <- file.path(dirs$main, "Run_nmf_profile.rds")
  writeLines("x", bundle)

  path <- NMFprofileR:::write_output_manifest(dirs, extra_paths = bundle)
  expect_true(file.exists(path))

  man <- readr::read_tsv(path, show_col_types = FALSE)
  expect_named(man, c("file", "type"))

  # The manifest never lists itself.
  expect_false(any(basename(man$file) == "manifest.tsv"))
  # Paths are relative to the results directory (no drive letters / leading sep).
  expect_false(any(grepl("^([A-Za-z]:)?/", man$file)))

  # Each representative file is classified as expected.
  type_of <- function(base) man$type[basename(man$file) == base]
  expect_equal(type_of("Run_nmf_profile.rds"), "profile_bundle")
  expect_equal(type_of("NMF_Result_Object_Rank_k2.rds"), "core_rds")
  expect_equal(type_of("Basis_Genes_Rank_k2.tsv"), "basis_genes")
  expect_equal(type_of("Marker_Genes_Rank_k2.tsv"), "marker_genes")
  expect_equal(type_of("Enrichment_Rank_k2_Factor_1.tsv"), "enrichment")
  expect_equal(type_of("Sample_Assignments_Rank_k2.tsv"), "sample_assignments")
  expect_equal(type_of("02_Rank_Survey_Plot.pdf"), "plot")
  expect_equal(type_of("Consolidated_Summary.tsv"), "summary")
  expect_equal(type_of("Run_run_log.txt"), "log")
})