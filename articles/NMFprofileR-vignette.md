# Getting started with NMFprofileR

## Overview

**NMFprofileR** runs a full multi-rank consensus Non-negative Matrix
Factorization (NMF) over a transcriptomics matrix and layers g:Profiler
functional enrichment on top. For each rank *k* it fits consensus NMF,
assigns genes and samples to factors, runs enrichment, computes
silhouette widths, and (optionally) writes a tree of tables and figures
to disk.

The heavy steps depend on the `NMF` package and on live g:Profiler
network calls, so those chunks below are shown but not executed. The
lightweight steps run for real.

``` r

library(NMFprofileR)
data("example_expression_data")
dim(example_expression_data)
#> [1] 20501    53
example_expression_data[1:5, 1:3]
#>       TCGA-02-0055 TCGA-06-0130 TCGA-06-0878
#> A1BG       391.804     152.8930     155.7130
#> A1CF         0.000       0.0000       0.0000
#> A2BP1      137.351      90.2527      15.8988
#> A2LD1       84.014      77.1377     181.5540
#> A2M      42876.300   22147.4000   14418.0000
```

## The one-call pipeline

[`NMFprofileR()`](https://milesgithub.github.io/NMFprofileR/reference/NMFprofileR.md)
orchestrates the whole workflow. It writes results under `output_prefix`
**and** returns an `nmf_profile` object with everything in memory.

``` r

result <- NMFprofileR(
  expression_data = example_expression_data,
  nmf_rank        = 2:5,
  output_prefix   = file.path(tempdir(), "NMF_demo"),
  nmf_nrun        = 20,
  gprofiler_sources = c("GO:BP", "GO:MF", "GO:CC", "REAC", "TF")
)
result
#> <nmf_profile>
#>   Ranks fitted : 2, 3, 4, 5
#>   Samples      : 53
#>   Runtime      : 3.2 mins
#>   g:Profiler   : e111_eg58_p18_...
#>   Output dir   : .../NMF_demo_Results
```

The result keeps the fitted objects, gene assignments, and enrichment
tables so you never have to re-read files:

``` r

summary(result)              # consolidated per-factor metrics across ranks
result$fits[["3"]]           # fitted NMF object for k = 3
result$basis_genes[["3"]]    # every gene -> its dominant factor (argmax)
result$marker_genes[["3"]]   # factor-specific marker genes (specificity-scored)
result$sample_assignments[["3"]]
result$enrichment$combined[["3"]]
```

Pass `write_files = FALSE` to compute entirely in memory without
creating any files – useful for scripting, examples, and tests.

## Preprocessing

The first stage validates and filters the matrix. It requires
gene-symbol rownames, resolves duplicated symbols (`on_duplicate_genes`,
default `"collapse_max"`), drops all-zero genes, and selects
highly-variable genes.

``` r

m <- nmf_preprocess(
  example_expression_data,
  expression_threshold = 10,
  variance_quantile    = 0.9
)
#> ! Dropping 794 all-zero gene row(s).
#> ℹ Filtering genes by mean expression < 10...
#> ℹ Filtering genes by variance < 0.9 quantile (raw scale)...
dim(m)
#> [1] 1574   53
```

`variance_scale = "log"` ranks variability on `log2(x + 1)` instead of
the raw scale, which reduces the mean-variance confound of count data.
The matrix fed to NMF is always on the original non-negative scale.

## Composable stages

Every stage is exported and works in memory, so you can build a custom
workflow or inspect intermediate results.

``` r

fit     <- nmf_fit(m, rank = 3, nrun = 20)
genes   <- nmf_basis_genes(fit)          # argmax assignment (all genes)
markers <- nmf_marker_genes(fit)         # specificity-scored markers (subset)
samples <- nmf_sample_assignments(fit)   # sample -> factor + silhouette

gene_sets <- split(genes$Gene, genes$Factor)
enr <- nmf_enrichment(gene_sets, background_genes = rownames(m))
```

## Choosing the number of components

[`nmf_rank_diagnostics()`](https://milesgithub.github.io/NMFprofileR/reference/nmf_rank_diagnostics.md)
lays the rank-selection metrics side by side. In a real run you would
pass `result$rank_metrics`; here is a small illustrative table:

``` r

rank_metrics <- data.frame(
  rank       = 2:5,
  cophenetic = c(0.98, 0.95, 0.90, 0.86),
  dispersion = c(0.90, 0.82, 0.75, 0.70),
  silhouette = c(0.61, 0.55, 0.48, 0.44)
)
nmf_rank_diagnostics(rank_metrics)
#> # A tibble: 12 × 3
#>     rank metric     value
#>    <dbl> <chr>      <dbl>
#>  1     2 cophenetic  0.98
#>  2     3 cophenetic  0.95
#>  3     4 cophenetic  0.9 
#>  4     5 cophenetic  0.86
#>  5     2 dispersion  0.9 
#>  6     3 dispersion  0.82
#>  7     4 dispersion  0.75
#>  8     5 dispersion  0.7 
#>  9     2 silhouette  0.61
#> 10     3 silhouette  0.55
#> 11     4 silhouette  0.48
#> 12     5 silhouette  0.44
```

A common heuristic is to pick the largest rank before the cophenetic
correlation drops sharply, cross-checked against a silhouette peak.

## Reproducibility

Every run records a manifest (parameters, package version, runtime), the
exact g:Profiler database version behind the enrichment, and full
session information, and the UMAP embedding is seeded – so a result can
be traced to the inputs and data snapshot that produced it
(`result$provenance`).

## Installation note

NMFprofileR depends on the Bioconductor packages `ComplexHeatmap` and
`circlize`. When installing from GitHub, enable the Bioconductor
repositories first:

``` r

setRepositories(ind = c(1, 2))   # CRAN + Bioconductor
devtools::install_github("MilesGithub/NMFprofileR")
```
