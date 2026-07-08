# Tidy NMF rank-selection diagnostics

Lays the key rank-selection metrics side by side across the tested ranks
so the "authoritative" component count can be chosen (typically the
largest cophenetic correlation before it drops, or a silhouette peak).
Works on the \`rank_metrics\` produced by \[NMFprofileR()\] (via
\`NMF::nmfEstimateRank\`).

## Usage

``` r
nmf_rank_diagnostics(x, metrics = c("cophenetic", "dispersion", "silhouette"))
```

## Arguments

- x:

  An \`nmf_profile\` object from \[NMFprofileR()\], or a
  \`rank_metrics\` data frame directly.

- metrics:

  Character vector of metrics to extract; matched case-insensitively and
  by prefix (so \`"silhouette"\` matches \`silhouette.consensus\`).
  Metrics not present in the data are silently skipped.

## Value

A tibble in long form with columns \`rank\`, \`metric\`, and \`value\`.

## Examples

``` r
rm <- data.frame(rank = 2:4, cophenetic = c(0.99, 0.95, 0.88),
                 dispersion = c(0.9, 0.8, 0.7), silhouette = c(0.6, 0.5, 0.4))
nmf_rank_diagnostics(rm)
#> # A tibble: 9 × 3
#>    rank metric     value
#>   <dbl> <chr>      <dbl>
#> 1     2 cophenetic  0.99
#> 2     3 cophenetic  0.95
#> 3     4 cophenetic  0.88
#> 4     2 dispersion  0.9 
#> 5     3 dispersion  0.8 
#> 6     4 dispersion  0.7 
#> 7     2 silhouette  0.6 
#> 8     3 silhouette  0.5 
#> 9     4 silhouette  0.4 
```
