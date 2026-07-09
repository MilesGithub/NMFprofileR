# Perform a Full Multi-Rank NMF and Functional Enrichment Analysis

This function is the main entry point for the NMFprofileR workflow. It
orchestrates the entire analysis pipeline, from preprocessing the input
expression matrix to performing a multi-rank NMF analysis, extracting
basis genes, running functional enrichment, and generating a
comprehensive suite of visualizations and summary reports.

## Usage

``` r
NMFprofileR(
  expression_data,
  nmf_rank = 2:5,
  output_prefix,
  gene_list_filter_file = NULL,
  expression_threshold = 10,
  variance_quantile = 0.7,
  nmf_method = "brunet",
  nmf_nrun = 20,
  nmf_seed = 123456,
  gprofiler_organism = "hsapiens",
  gprofiler_sources = c("GO:BP", "GO:CC", "GO:MF", "REAC", "TF"),
  gprofiler_correction = "g_SCS",
  gprofiler_cutoff = 0.05,
  enrichment_plot_top_n = 50,
  verbose = FALSE,
  write_files = TRUE,
  on_duplicate_genes = c("collapse_max", "collapse_mean", "error"),
  emit_marker_genes = FALSE,
  variance_scale = c("raw", "log"),
  custom_theme = NULL,
  factor_palette = NULL,
  umap_n_neighbors = 15,
  run_id = NULL,
  nmf_parallel = FALSE,
  nmf_cores = NULL,
  max_query_size = 10000,
  enrichment_cache = NULL
)
```

## Arguments

- expression_data:

  A data frame or matrix of normalized, non-negative gene expression
  values. Rows should represent genes (with gene symbols as rownames)
  and columns should represent samples.

- nmf_rank:

  An integer vector specifying the factorization ranks (k) to test
  (e.g., \`2:5\`).

- output_prefix:

  A character string defining the base path and file name prefix for all
  output files and directories. For example, \`"results/MyAnalysis"\`.

- gene_list_filter_file:

  An optional character string path to a plain text file containing a
  list of gene symbols (one per line). If provided, the expression
  matrix will be pre-filtered to include only these genes.

- expression_threshold:

  A numeric value. Genes with a mean expression across all samples below
  this threshold will be removed prior to NMF.

- variance_quantile:

  A numeric value between 0 and 1. After mean-based filtering, genes
  with variance below this quantile are removed. This step enriches for
  more informative genes.

- nmf_method:

  A character string specifying the NMF algorithm to use. Defaults to
  'brunet'. See the \`NMF::nmf\` documentation for other options.

- nmf_nrun:

  An integer specifying the number of runs for the consensus NMF. A
  higher number (e.g., 50-100) provides more stable results.

- nmf_seed:

  An integer to use as the random seed for reproducibility.

- gprofiler_organism:

  A character string specifying the organism for g:Profiler enrichment
  analysis (e.g., 'hsapiens', 'mmusculus').

- gprofiler_sources:

  A character vector of data sources to query in g:Profiler (e.g.,
  c("GO:BP", "REAC", "KEGG")).

- gprofiler_correction:

  A character string specifying the multiple-testing correction method
  used by g:Profiler. Defaults to \`"g_SCS"\`, g:Profiler's native Set
  Counts and Sizes method, which accounts for the overlapping,
  hierarchical structure of GO and pathway terms and is generally more
  appropriate (and less conservative) there than Bonferroni or
  Benjamini- Hochberg FDR. The query uses a custom background of the
  genes that survive preprocessing (\`custom_bg\`), so enrichment is
  assessed against the tested gene universe rather than the whole
  genome. Alternatives (\`"fdr"\`, \`"bonferroni"\`) are accepted for
  reviewer familiarity at some loss of power.

- gprofiler_cutoff:

  A numeric value for the significance threshold for enrichment results.

- enrichment_plot_top_n:

  An integer specifying the number of top enriched terms to display in
  the summary dot plots.

- verbose:

  A logical value. If TRUE, prints detailed messages. Defaults to FALSE.

- write_files:

  A logical value. If TRUE (the default), all results, plots, and log
  files are written to disk under \`output_prefix\`. If FALSE, nothing
  is written and no directories are created; the pipeline computes in
  memory and returns its results only. Useful for programmatic use,
  examples, and tests (which must not write outside \`tempdir()\`).

- on_duplicate_genes:

  How to handle duplicated gene symbols in the rownames of
  \`expression_data\`: \`"collapse_max"\` (the default) keeps the
  per-sample maximum across duplicate rows (the most expressed probe),
  \`"collapse_mean"\` averages them, and \`"error"\` stops with a
  message.

- emit_marker_genes:

  A logical value. If TRUE, in addition to the argmax basis-gene
  assignment the pipeline also extracts factor-specific marker genes
  (\`NMF::extractFeatures\`, Kim-Park specificity) and runs enrichment
  on them, emitting \`Marker_Genes\_\*\` and \`Marker_Enrichment\_\*\`
  outputs alongside the argmax results. Defaults to FALSE: markers are
  opt-in because computing them runs a second g:Profiler query per rank
  (roughly doubling the network work), which matters most for batch use.
  The argmax assignment is always the primary output. Regardless of this
  setting the return value always carries the \`marker_genes\` and
  \`enrichment\$markers\` elements (empty when FALSE).

- variance_scale:

  Either \`"raw"\` (default) or \`"log"\`. The scale on which gene
  variance is measured for the highly-variable-gene filter.
  \`expression_data\` is expected to be normalized, non-negative values
  (e.g. normalized counts or TPM); on such data \`expression_threshold =
  10\` and raw-scale variance are meaningful. \`"log"\` ranks
  variability on \`log2(x + 1)\`, reducing the mean-variance confound of
  count data. The matrix fed to NMF is always the original (raw,
  non-negative) scale regardless of this setting.

- custom_theme:

  An optional \`ggplot2\` theme object applied to the ggplot-based plots
  (the factor-summary bar plot, enrichment dot plots, and the sample
  UMAP) in place of the package's built-in theme. \`NULL\` (the default)
  uses the built-in theme.

- factor_palette:

  An optional character vector of colours used to colour NMF factors
  across the plots. \`NULL\` (the default) uses the built-in palette.
  When a rank has more factors than supplied colours, the palette is
  extended with \`grDevices::hcl.colors()\`.

- umap_n_neighbors:

  Integer; the number of neighbours for the sample-coefficient UMAP
  embedding (default 15). The UMAP is only drawn when the cohort has
  more samples than this value, and the effective neighbour count is
  capped at \`n_samples - 1\` so smaller cohorts still embed.

- run_id:

  An optional identifier for this run. When supplied it is stamped as a
  leading \`Run_ID\` column in \`consolidated_summary_df\` (and the
  on-disk \`Consolidated_Summary.tsv\`) and recorded in the run
  manifest, so results from many runs can be pooled and traced back to
  their source.

- nmf_parallel:

  A logical value. If \`FALSE\` (the default) the per-rank consensus
  runs are fitted sequentially. If \`TRUE\` they are spread across a
  local cluster; see \[nmf_fit()\] for details. For a fixed \`nmf_seed\`
  the parallel path yields the same factorization as the sequential one.
  The rank-estimation survey always runs sequentially.

- nmf_cores:

  Number of worker processes to use when \`nmf_parallel = TRUE\`.
  \`NULL\` (the default) uses one fewer than the number of detected
  cores.

- max_query_size:

  An integer. Factor gene sets larger than this are skipped (with a
  warning) rather than sent to g:Profiler, guarding against pathological
  oversized queries. Defaults to 10000.

- enrichment_cache:

  An optional path to a directory used to cache g:Profiler results. When
  set, each query is keyed by a hash of its genes, background, sources,
  and settings and read from / written to this directory, so re-runs
  (e.g. across a batch) skip the network. \`NULL\` (the default)
  disables caching.

## Value

Invisibly, an S3 object of class \`nmf_profile\` (a list underneath, so
\`\$\` access works) with \`print()\` and \`summary()\` methods.
Elements:

- \`consolidated_summary_df\`:

  A data frame summarizing every factor from every rank tested.

- \`rank_metrics\`:

  A data frame of quality metrics from the \`NMF::nmfEstimateRank\` rank
  survey (empty if the survey failed).

- \`fits\`:

  A named list of the fitted \`NMFfit\` objects, keyed by rank.

- \`basis_genes\`:

  A named list of per-rank basis-gene tables (every gene to its dominant
  factor, by argmax).

- \`marker_genes\`:

  A named list of per-rank factor-specific marker-gene tables
  (specificity-scored; only when \`emit_marker_genes = TRUE\`).

- \`sample_assignments\`:

  A named list of per-rank sample-assignment tables (sample to dominant
  factor, with silhouette widths).

- \`enrichment\`:

  A list with \`per_factor\`, \`combined\`, and \`markers\` named lists
  of per-rank g:Profiler enrichment data frames.

- \`nmf_rds\`:

  A named list of file paths to the saved NMF fit objects (empty when
  \`write_files = FALSE\`).

- \`output_dirs\`:

  A named list of the output directories created for the run, or
  \`NULL\` when \`write_files = FALSE\`.

- \`failures\`:

  A data frame (columns \`Rank\`, \`Reason\`) of ranks whose NMF fit
  failed and were skipped; empty when every rank succeeded.

- \`runtime\`:

  A \`difftime\` giving the total pipeline runtime.

- \`provenance\`:

  A named list recording the captured g:Profiler database version
  (\`gprofiler_version\`), the \`run_parameters\` used, and the paths of
  the manifest and session-information files written (\`files\`), for
  reproducibility.

Unless \`write_files = FALSE\`, all detailed results, plots, and log
files are additionally written to disk under the location specified by
\`output_prefix\`.

## Reproducibility

NMF is run sequentially (not across a parallel backend), which is both
what keeps results reproducible for a fixed \`nmf_seed\` and what avoids
a failure mode of NMF's parallel execution when called from within a
package. For a fixed \`nmf_seed\` the factorization is therefore
deterministic within a fixed computational environment; exact numerical
reproducibility across machines can still be affected by the BLAS
implementation and threading. Every run also writes a manifest and
session information (surfaced in the \`provenance\` element of the
return value) so a result can be traced to the parameters and g:Profiler
database snapshot that produced it.
