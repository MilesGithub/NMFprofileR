# Extract factor-specific marker genes by specificity score

Unlike \[nmf_basis_genes()\] (which assigns \*every\* gene to its
dominant factor by argmax), this selects only the genes that are
specifically associated with each factor, using
\`NMF::extractFeatures()\` (the Kim & Park specificity scoring). The
result is a smaller, sharper set of marker genes per factor; many genes
are intentionally left unassigned.

## Usage

``` r
nmf_marker_genes(fit)
```

## Arguments

- fit:

  An \`NMFfit\` object from \[nmf_fit()\].

## Value

A tibble with columns \`Gene\` and \`Factor\`. Factors with no specific
markers contribute no rows.

## See also

\[nmf_basis_genes()\], \[nmf_enrichment()\]

## Examples

``` r
if (FALSE) { # \dontrun{
markers <- nmf_marker_genes(fit)
split(markers$Gene, markers$Factor)
} # }
```
