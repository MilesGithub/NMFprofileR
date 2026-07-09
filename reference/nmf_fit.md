# Fit consensus NMF at a single rank

A thin, in-memory wrapper around \[NMF::nmf()\] with the package's
platform-aware parallel handling and error trapping. Returns the fit
object without writing anything to disk.

## Usage

``` r
nmf_fit(
  expr_matrix,
  rank,
  method = "brunet",
  nrun = 20,
  seed = 123456,
  nmf_parallel = FALSE,
  n_cores = NULL
)
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

- nmf_parallel:

  Logical; if \`FALSE\` (the default) the \`nrun\` consensus runs are
  executed sequentially. If \`TRUE\` they are spread across a local
  PSOCK cluster. The parallel path is opt-in because NMF's default
  parallel backend fails when NMF is invoked from inside another package
  (its workers cannot see the NMF namespace); here the cluster is
  configured to load NMF on each worker, which both fixes that and,
  because NMF pre-generates the RNG stream for every run, yields results
  identical to the sequential path for a fixed \`seed\`. Falls back to
  sequential (with a warning) if the cluster cannot be created.

- n_cores:

  Number of worker processes to use when \`nmf_parallel = TRUE\`.
  \`NULL\` (the default) uses one fewer than the number of detected
  cores, capped at \`nrun\`.

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
