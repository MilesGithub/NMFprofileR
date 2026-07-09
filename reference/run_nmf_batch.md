# Run NMFprofileR over many cohorts and consolidate the results

A batch driver that applies \[NMFprofileR()\] to each of several
expression matrices ("cohorts"), writing each cohort's outputs to its
own subdirectory of \`output_dir\` and returning a single consolidated
per-factor summary across all cohorts. It is resilient (one cohort
failing does not abort the batch), resumable (already-completed cohorts
can be skipped), and tags every row with the cohort's \`run_id\` so
results can be pooled and traced.

## Usage

``` r
run_nmf_batch(
  cohorts,
  output_dir,
  ...,
  skip_existing = TRUE,
  on_error = c("continue", "stop"),
  write_batch_summary = TRUE
)
```

## Arguments

- cohorts:

  A named list of expression matrices or data frames (one per cohort),
  in the form accepted by \[NMFprofileR()\]. The names are used both as
  the \`run_id\` and as the output subdirectory for each cohort.

- output_dir:

  A directory under which each cohort's \`\<name\>\_Results\` tree is
  written.

- ...:

  Further arguments passed to \[NMFprofileR()\] (e.g. \`nmf_rank\`,
  \`nmf_nrun\`, \`gprofiler_sources\`, \`nmf_parallel\`,
  \`enrichment_cache\`). Do not pass \`expression_data\`,
  \`output_prefix\`, or \`run_id\`; the driver sets those per cohort.

- skip_existing:

  Logical. If \`TRUE\` (the default), a cohort whose results directory
  already exists is not re-run; its summary is loaded from disk and
  still included in the consolidated output.

- on_error:

  One of \`"continue"\` (the default; log the failing cohort and carry
  on) or \`"stop"\` (abort the whole batch on the first cohort error).

- write_batch_summary:

  Logical. If \`TRUE\` (the default) the consolidated summary is also
  written to \`output_dir/Batch_Consolidated_Summary.tsv\`.

## Value

Invisibly, a list with:

- \`consolidated\`:

  A data frame binding every cohort's \`consolidated_summary_df\` (each
  tagged with its \`Run_ID\`).

- \`failures\`:

  A data frame (\`Cohort\`, \`Reason\`) of cohorts that errored.

- \`skipped\`:

  A character vector of cohorts skipped because their results already
  existed.

## See also

\[NMFprofileR()\]

## Examples

``` r
if (FALSE) { # \dontrun{
data("example_expression_data")
# Split the samples into two illustrative cohorts.
half <- ncol(example_expression_data) %/% 2
cohorts <- list(
  cohort_A = example_expression_data[, seq_len(half)],
  cohort_B = example_expression_data[, (half + 1):ncol(example_expression_data)]
)
batch <- run_nmf_batch(cohorts, output_dir = "batch_out", nmf_rank = 2:3, nmf_nrun = 10)
batch$consolidated
} # }
```
