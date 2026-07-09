# Batch analysis across many cohorts

## Overview

A common workflow is to run the same NMF profiling on many independent
expression matrices – different cohorts, subgroups, tissues, or
conditions – and then compare the factors found across them.
[`run_nmf_batch()`](https://milesgithub.github.io/NMFprofileR/reference/run_nmf_batch.md)
drives
[`NMFprofileR()`](https://milesgithub.github.io/NMFprofileR/reference/NMFprofileR.md)
over a named list of cohorts and returns a single consolidated summary
tagged by cohort, so the results can be pooled and traced back to their
source.

The driver is:

- **resilient** – one cohort failing does not abort the batch
  (`on_error = "continue"`);
- **resumable** – cohorts whose results already exist are skipped and
  loaded from disk (`skip_existing = TRUE`);
- **self-describing** – every summary row carries the cohort’s `Run_ID`.

``` r

library(NMFprofileR)
data("example_expression_data")
dim(example_expression_data)
#> [1] 20501    53
```

## Defining cohorts

A cohort is just an expression matrix in the form
[`NMFprofileR()`](https://milesgithub.github.io/NMFprofileR/reference/NMFprofileR.md)
accepts (genes in rows with gene-symbol rownames, samples in columns).
`cohorts` is a **named** list; the names become both the `run_id` and
the output subdirectory for each cohort.

Here we split the bundled example data into two illustrative cohorts:

``` r

half <- ncol(example_expression_data) %/% 2
cohorts <- list(
  cohort_A = example_expression_data[, seq_len(half)],
  cohort_B = example_expression_data[, (half + 1):ncol(example_expression_data)]
)
vapply(cohorts, ncol, integer(1))
#> cohort_A cohort_B 
#>       26       27
```

In practice the cohorts might come from splitting one matrix by a
grouping variable, or from reading several matrices off disk into a
named list.

## Running the batch

Every argument other than the ones the driver manages
(`expression_data`, `output_prefix`, `run_id`) is passed straight
through to
[`NMFprofileR()`](https://milesgithub.github.io/NMFprofileR/reference/NMFprofileR.md).
Each cohort’s outputs land in `output_dir/<name>_Results/`, and the
consolidated summary is also written to
`output_dir/Batch_Consolidated_Summary.tsv`.

``` r

batch <- run_nmf_batch(
  cohorts,
  output_dir        = file.path(tempdir(), "batch_demo"),
  nmf_rank          = 2:4,
  nmf_nrun          = 20,
  gprofiler_sources = c("GO:BP", "REAC"),
  enrichment_cache  = file.path(tempdir(), "gprofiler_cache")
)

batch$consolidated   # one row per factor, per rank, per cohort (tagged Run_ID)
batch$skipped        # cohorts skipped because results already existed
batch$failures       # cohorts that errored (Cohort, Reason)
```

The consolidated table has the same per-factor columns as a single run,
with a leading `Run_ID`:

    #>     Run_ID Rank   Factor Num_Samples Percent_Samples
    #> 1 cohort_A    2 Factor_1          14            53.8
    #> 2 cohort_A    2 Factor_2          12            46.2
    #> 3 cohort_B    2 Factor_1          15            57.7

## Resuming an interrupted batch

Because `skip_existing = TRUE` by default, re-running the same call
after an interruption only computes the cohorts that have not finished
yet; the completed ones are read back from their saved `nmf_profile`
bundle (or `Consolidated_Summary.tsv`) and still included in
`batch$consolidated`. Set `skip_existing = FALSE` to force every cohort
to re-run.

## Making a batch fast and robust

- **Cache enrichment.** Pass `enrichment_cache = <dir>` so repeated
  g:Profiler queries (common when cohorts share genes, or when
  re-running) are served from disk instead of the network. Queries are
  keyed by a hash of their genes, background, sources, and settings.
- **Parallelize the fits.** Pass `nmf_parallel = TRUE` (optionally with
  `nmf_cores`) to spread each cohort’s consensus runs across a local
  cluster. For a fixed `nmf_seed` the parallel path returns the same
  factorization as the sequential one.
- **Keep going on failure.** With `on_error = "continue"` (the default)
  a cohort that errors is recorded in `batch$failures` and the batch
  proceeds; use `on_error = "stop"` to fail fast.
- **Skip the second enrichment pass.** Marker-gene enrichment is opt-in
  (`emit_marker_genes = FALSE` by default), which halves the number of
  g:Profiler queries per cohort.

## Notes

The heavy steps depend on the `NMF` package and on live g:Profiler
network calls, so the batch chunk above is shown but not executed; the
cohort-building chunks run for real.

``` r

utils::sessionInfo()
```
