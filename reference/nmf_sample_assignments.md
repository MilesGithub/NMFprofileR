# Assign samples to their dominant NMF factor with silhouette widths

Assigns every sample to its dominant coefficient factor (argmax over the
coefficient matrix) and attaches the consensus-matrix silhouette width
for each sample.

## Usage

``` r
nmf_sample_assignments(fit, verbose = FALSE)
```

## Arguments

- fit:

  An \`NMFfit\` object from \[nmf_fit()\].

- verbose:

  Logical; passed to the silhouette computation for diagnostic messages.

## Value

A data frame with columns \`SampleID\`, \`Dominant_Factor\`, and
\`Silhouette_NMF\`, ordered by factor then descending silhouette.

## See also

\[nmf_fit()\]

## Examples

``` r
if (FALSE) { # \dontrun{
nmf_sample_assignments(fit)
} # }
```
