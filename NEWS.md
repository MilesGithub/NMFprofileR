# NMFprofileR 0.4.0

## New features

* `run_nmf_batch()` runs `NMFprofileR()` over a named list of cohorts, writing
  each to its own subdirectory and returning a single consolidated per-factor
  summary tagged by cohort `Run_ID`. It is resilient (`on_error = "continue"`),
  resumable (`skip_existing` reloads finished cohorts from disk), and writes a
  `Batch_Consolidated_Summary.tsv`.
* Opt-in parallel NMF: `nmf_parallel = TRUE` (with optional `nmf_cores`) spreads
  the consensus runs across a local cluster. The cluster loads NMF on each
  worker, so it works when called from inside a package, and for a fixed
  `nmf_seed` it returns the same factorization as the sequential path. The
  default remains sequential.
* g:Profiler hardening: an oversized-query guard (`max_query_size`), automatic
  retry with backoff on transient network failures, and an optional on-disk
  result cache (`enrichment_cache`) keyed by query hash so re-runs and batches
  skip the network.

## Infrastructure

* New "Batch analysis across many cohorts" vignette.

# NMFprofileR 0.3.0

## Breaking changes

* `emit_marker_genes` now defaults to `FALSE` (previously `TRUE`). Specificity
  marker genes and their enrichment are now opt-in, halving the g:Profiler
  network work of a default run. Set `emit_marker_genes = TRUE` to restore the
  previous behaviour. The return value still always carries the `marker_genes`
  and `enrichment$markers` elements (empty when markers are off).

## New features

* Themeable plots: `NMFprofileR()` gains `custom_theme` (a `ggplot2` theme
  applied to the ggplot-based plots) and `factor_palette` (a colour vector for
  factors), so figures can match a project's house style without editing the
  package.
* `umap_n_neighbors` parametrizes the sample-coefficient UMAP. The embedding is
  drawn whenever the cohort has more samples than this value, and the neighbour
  count is capped at `n_samples - 1`, so smaller cohorts now get a UMAP too.
* Self-describing outputs: an optional `run_id` is stamped as a leading `Run_ID`
  column in the consolidated summary and recorded in the manifest; the whole
  `nmf_profile` result is saved as a single `<prefix>_nmf_profile.rds` bundle;
  and a `Summaries/manifest.tsv` lists every file a run produced with its type.
* Failed ranks are reported: the result gains a `failures` data frame (rank and
  reason) and `print()` lists any failed ranks, instead of ranks silently
  vanishing from the output.

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