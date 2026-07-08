# Fit consensus NMF at a single rank

A thin, in-memory wrapper around \[NMF::nmf()\] with the package's
platform-aware parallel handling and error trapping. Returns the fit
object without writing anything to disk.

## Usage

``` r
nmf_fit(expr_matrix, rank, method = "brunet", nrun = 20, seed = 123456)
```

## Arguments

- expr_matrix:

  A preprocessed, non-negative numeric matrix (see
  \[nmf_preprocess()\]).

- rank:

  The factorization rank (integer).

- method:

  The NMF algorithm (e.g. \`"brunet"\`).

- nrun:

  The number of consensus runs.

- seed:

  An integer seed for reproducibility.

## Value

An \`NMFfit\` object, or \`NULL\` if the fit fails.

## See also

\[nmf_preprocess()\], \[nmf_basis_genes()\],
\[nmf_sample_assignments()\]

## Examples

``` r
if (FALSE) { # \dontrun{
m <- nmf_preprocess(example_expression_data, expression_threshold = 0,
                    variance_quantile = 0)
fit <- nmf_fit(m, rank = 3, nrun = 10)
} # }
```
