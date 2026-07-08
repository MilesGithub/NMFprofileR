# NMFprofileR 0.2.0

## New features

* The pipeline is now available as composable, exported stages that return
  results in memory: `nmf_preprocess()`, `nmf_fit()`, `nmf_basis_genes()`,
  `nmf_marker_genes()`, `nmf_sample_assignments()`, `nmf_enrichment()`, and
  `nmf_rank_diagnostics()`. `NMFprofileR()` orchestrates them.
* `NMFprofileR()` returns an S3 `nmf_profile` object (with `print()` and
  `summary()` methods) carrying the fitted objects, gene assignments, and
  enrichment tables in memory, not just file paths.
* `write_files = FALSE` runs the pipeline entirely in memory, creating no files.
* Up-front input validation with actionable errors, plus explicit duplicate
  gene-symbol handling via `on_duplicate_genes` (default `"collapse_max"`).
* Factor-specific marker genes (`NMF::extractFeatures`, Kim-Park specificity)
  are emitted alongside the argmax basis-gene assignment (`emit_marker_genes`,
  default `TRUE`).
* `variance_scale` selects raw- or log-scale variance for highly-variable-gene
  selection (default `"raw"`, preserving previous behaviour).
* Reproducibility hardening: each run records a manifest, the g:Profiler
  database version, and session information; the UMAP embedding is seeded.

## Bug fixes and consistency

* The global expression heatmap now groups genes by the same assignment used for
  the exported basis-gene table (previously the figure and table could disagree).
* Cleared ggplot2 (>= 3.4) `linewidth` deprecation warnings.
* `capture_plot()` no longer leaks a graphics device if a plot expression fails.
* Clear error (instead of an opaque NMF crash) when fewer than two genes survive
  filtering.

## Infrastructure

* Added a testthat (edition 3) test suite, modernized the R-CMD-check workflow
  (r-lib/actions matrix), and added a test-coverage workflow.

# NMFprofileR 0.1.0

* Initial version.