# Assess NMF factor stability by subsampling

Quantifies how reproducible each factor is by refitting NMF on random
subsets of the samples and measuring how well the resampled factors
recover the reference factors (fitted on the full matrix). For each
resample the basis loadings are correlated against the reference and
each reference factor's best match is recorded; the per-factor stability
is the mean of these best matches across resamples (1 = perfectly
reproducible).

## Usage

``` r
nmf_stability(
  expr_matrix,
  rank,
  method = "brunet",
  nrun = 10,
  seed = 123456,
  nboot = 30,
  subsample_frac = 0.8,
  nmf_parallel = FALSE,
  n_cores = NULL
)
```

## Arguments

- expr_matrix:

  A preprocessed, non-negative numeric matrix (genes x samples), as from
  \[nmf_preprocess()\].

- rank:

  The factorization rank (integer).

- method, nrun, seed:

  Passed to \[nmf_fit()\] for every fit.

- nboot:

  Number of subsampling iterations (default 30).

- subsample_frac:

  Fraction of samples drawn (without replacement) in each iteration
  (default 0.8).

- nmf_parallel, n_cores:

  Passed to \[nmf_fit()\].

## Value

A tibble with columns \`Factor\` and \`Stability\` (mean best-match
correlation across resamples), one row per reference factor.

## See also

\[nmf_fit()\], \[nmf_rank_diagnostics()\]

## Examples

``` r
if (FALSE) { # \dontrun{
m <- nmf_preprocess(example_expression_data, expression_threshold = 10,
                    variance_quantile = 0.9)
nmf_stability(m, rank = 3, nrun = 10, nboot = 20)
} # }
```
