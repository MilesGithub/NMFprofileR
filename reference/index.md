# Package index

## Pipeline

The end-to-end orchestrator and batch driver.

- [`NMFprofileR()`](https://milesgithub.github.io/NMFprofileR/reference/NMFprofileR.md)
  : Perform a Full Multi-Rank NMF and Functional Enrichment Analysis
- [`run_nmf_batch()`](https://milesgithub.github.io/NMFprofileR/reference/run_nmf_batch.md)
  : Run NMFprofileR over many cohorts and consolidate the results

## Composable stages

Individual, in-memory pipeline steps.

- [`nmf_preprocess()`](https://milesgithub.github.io/NMFprofileR/reference/nmf_preprocess.md)
  : Preprocess an expression matrix for NMF
- [`nmf_fit()`](https://milesgithub.github.io/NMFprofileR/reference/nmf_fit.md)
  : Fit consensus NMF at a single rank
- [`nmf_basis_genes()`](https://milesgithub.github.io/NMFprofileR/reference/nmf_basis_genes.md)
  : Assign genes to their dominant NMF factor
- [`nmf_marker_genes()`](https://milesgithub.github.io/NMFprofileR/reference/nmf_marker_genes.md)
  : Extract factor-specific marker genes by specificity score
- [`nmf_sample_assignments()`](https://milesgithub.github.io/NMFprofileR/reference/nmf_sample_assignments.md)
  : Assign samples to their dominant NMF factor with silhouette widths
- [`nmf_enrichment()`](https://milesgithub.github.io/NMFprofileR/reference/nmf_enrichment.md)
  : Run g:Profiler enrichment on NMF factor gene sets
- [`nmf_rank_diagnostics()`](https://milesgithub.github.io/NMFprofileR/reference/nmf_rank_diagnostics.md)
  : Tidy NMF rank-selection diagnostics

## Result object

The nmf_profile S3 object.

- [`print(`*`<nmf_profile>`*`)`](https://milesgithub.github.io/NMFprofileR/reference/print.nmf_profile.md)
  : Print an NMFprofileR result
- [`summary(`*`<nmf_profile>`*`)`](https://milesgithub.github.io/NMFprofileR/reference/summary.nmf_profile.md)
  : Summarize an NMFprofileR result

## Data

- [`example_expression_data`](https://milesgithub.github.io/NMFprofileR/reference/example_expression_data.md)
  : Example Expression Data for NMFprofileR
