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
  emit_marker_genes = TRUE,
  variance_scale = c("raw", "log")
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

  A logical value. If TRUE (the default), in addition to the argmax
  basis-gene assignment the pipeline also extracts factor-specific
  marker genes (\`NMF::extractFeatures\`, Kim-Park specificity) and runs
  enrichment on them, emitting \`Marker_Genes\_\*\` and
  \`Marker_Enrichment\_\*\` outputs alongside the argmax results. The
  argmax assignment remains the primary, unchanged output.

- variance_scale:

  Either \`"raw"\` (default) or \`"log"\`. The scale on which gene
  variance is measured for the highly-variable-gene filter.
  \`expression_data\` is expected to be normalized, non-negative values
  (e.g. normalized counts or TPM); on such data \`expression_threshold =
  10\` and raw-scale variance are meaningful. \`"log"\` ranks
  variability on \`log2(x + 1)\`, reducing the mean-variance confound of
  count data. The matrix fed to NMF is always the original (raw,
  non-negative) scale regardless of this setting.

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

For a fixed \`nmf_seed\` the NMF factorization is deterministic and
yields identical basis and coefficient matrices whether NMF runs
sequentially or across a parallel backend, within a fixed computational
environment. Exact numerical reproducibility across machines can still
be affected by the BLAS implementation and threading. Every run also
writes a manifest and session information (surfaced in the
\`provenance\` element of the return value) so a result can be traced to
the parameters and g:Profiler database snapshot that produced it.
