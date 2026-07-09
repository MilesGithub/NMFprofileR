# Align NMF factors across runs by basis-loading correlation

Matches the factors of several NMF runs to a common reference so that
"the same" factor can be tracked across runs (e.g. across cohorts,
ranks, or resamples). The first element of \`x\` is the reference; every
other run's factors are correlated against the reference factors on
their shared genes and matched one-to-one.

## Usage

``` r
align_factors(
  x,
  method = c("greedy", "hungarian"),
  min_cor = 0.5,
  cor_method = c("pearson", "spearman")
)
```

## Arguments

- x:

  A named list (length \>= 2) of \`NMFfit\` objects and/or basis
  matrices (genes in rows with gene-symbol rownames, factors in
  columns). The names identify the runs; the first element is the
  reference.

- method:

  How to pick the one-to-one matches: \`"greedy"\` (the default)
  repeatedly takes the highest remaining correlation; \`"hungarian"\`
  finds the globally optimal assignment via \`clue::solve_LSAP()\` (used
  only if the \`clue\` package is installed, otherwise it falls back to
  greedy).

- min_cor:

  Minimum absolute correlation for a pair to be reported as a match
  (default 0.5).

- cor_method:

  Correlation method, \`"pearson"\` (default) or \`"spearman"\`.

## Value

A list with:

- \`reference\`:

  The name of the reference run.

- \`correlations\`:

  A tibble (\`run\`, \`factor\`, \`ref_factor\`, \`correlation\`) of
  every run factor against every reference factor.

- \`matches\`:

  A tibble (\`run\`, \`factor\`, \`ref_factor\`, \`correlation\`) of the
  selected one-to-one matches with \`abs(correlation) \>= min_cor\`.

## See also

\[nmf_fit()\]

## Examples

``` r
if (FALSE) { # \dontrun{
fit_a <- nmf_fit(m, rank = 3)
fit_b <- nmf_fit(m, rank = 3, seed = 99)
aln <- align_factors(list(run_a = fit_a, run_b = fit_b))
aln$matches
} # }
```
