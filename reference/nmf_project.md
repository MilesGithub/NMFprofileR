# Project new samples onto learned NMF factors

Scores new samples against an existing factorization by solving for
their coefficient matrix H with the basis matrix W held fixed (\`min
\|\|V - W H\|\|\` s.t. \`H \>= 0\`), using non-negative multiplicative
updates. Genes are aligned on the shared symbols between W and the new
data.

## Usage

``` r
nmf_project(fit, newdata, n_iter = 200L, tol = 1e-06)
```

## Arguments

- fit:

  An \`NMFfit\` object (or a numeric basis matrix W, genes x factors)
  defining the learned factors.

- newdata:

  A numeric matrix or data frame of new samples (genes in rows with
  gene-symbol rownames, samples in columns).

- n_iter:

  Maximum number of multiplicative-update iterations (default 200).

- tol:

  Convergence tolerance on the maximum change in H (default 1e-6).

## Value

A data frame with one row per new sample: \`SampleID\`, one coefficient
column per factor (\`Factor_1 ... Factor_k\`), and \`Dominant_Factor\`.

## See also

\[nmf_fit()\], \[nmf_sample_assignments()\]

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- nmf_fit(m, rank = 3)
nmf_project(fit, newdata = m_new)
} # }
```
