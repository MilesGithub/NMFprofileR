# Assign genes to their dominant NMF factor

Assigns every gene to its dominant basis factor (argmax over the basis
matrix, via \[NMF::predict()\]), returning a tidy table ordered by
factor then gene.

## Usage

``` r
nmf_basis_genes(fit)
```

## Arguments

- fit:

  An \`NMFfit\` object from \[nmf_fit()\].

## Value

A tibble with columns \`Gene\` and \`Factor\` (a factor with levels
\`Factor_1 ... Factor_k\`).

## See also

\[nmf_fit()\], \[nmf_enrichment()\]

## Examples

``` r
if (FALSE) { # \dontrun{
basis_genes <- nmf_basis_genes(fit)
split(basis_genes$Gene, basis_genes$Factor)
} # }
```
