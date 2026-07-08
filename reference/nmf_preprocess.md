# Preprocess an expression matrix for NMF

Validates the input (numeric matrix, gene-symbol rownames,
duplicate-symbol resolution, all-zero row/column checks) and applies the
standard filtering (optional gene-list filter, mean-expression
threshold, variance quantile, non-negativity). This is the first stage
of the \[NMFprofileR()\] pipeline, exposed for composable use.

## Usage

``` r
nmf_preprocess(
  expression_data,
  expression_threshold = 10,
  variance_quantile = 0.7,
  gene_list_filter_file = NULL,
  on_duplicate_genes = c("collapse_max", "collapse_mean", "error"),
  variance_scale = c("raw", "log")
)
```

## Arguments

- expression_data:

  A data frame or matrix of normalized, non-negative gene expression
  values. Rows should represent genes (with gene symbols as rownames)
  and columns should represent samples.

- expression_threshold:

  A numeric value. Genes with a mean expression across all samples below
  this threshold will be removed prior to NMF.

- variance_quantile:

  A numeric value between 0 and 1. After mean-based filtering, genes
  with variance below this quantile are removed. This step enriches for
  more informative genes.

- gene_list_filter_file:

  An optional character string path to a plain text file containing a
  list of gene symbols (one per line). If provided, the expression
  matrix will be pre-filtered to include only these genes.

- on_duplicate_genes:

  How to handle duplicated gene symbols in the rownames of
  \`expression_data\`: \`"collapse_max"\` (the default) keeps the
  per-sample maximum across duplicate rows (the most expressed probe),
  \`"collapse_mean"\` averages them, and \`"error"\` stops with a
  message.

- variance_scale:

  Either \`"raw"\` (default) or \`"log"\`. The scale on which gene
  variance is measured for the highly-variable-gene filter.
  \`expression_data\` is expected to be normalized, non-negative values
  (e.g. normalized counts or TPM); on such data \`expression_threshold =
  10\` and raw-scale variance are meaningful. \`"log"\` ranks
  variability on \`log2(x + 1)\`, reducing the mean-variance confound of
  count data. The matrix fed to NMF is always the original (raw,
  non-negative) scale regardless of this setting.

## Value

A filtered, non-negative numeric matrix ready for \[nmf_fit()\].

## See also

\[nmf_fit()\], \[NMFprofileR()\]

## Examples

``` r
data("example_expression_data")
m <- nmf_preprocess(example_expression_data, expression_threshold = 0,
                    variance_quantile = 0)
#> ! Dropping 794 all-zero gene row(s).
#> ℹ Filtering genes by mean expression < 0...
#> ℹ Filtering genes by variance < 0 quantile (raw scale)...
dim(m)
#> [1] 19707    53
```
