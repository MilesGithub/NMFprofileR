# Run g:Profiler enrichment on NMF factor gene sets

Runs per-factor and combined g:Profiler functional enrichment on a list
of factor gene sets, in memory (nothing is written to disk).

## Usage

``` r
nmf_enrichment(
  basis_genes_list,
  background_genes,
  organism = "hsapiens",
  sources = c("GO:BP", "GO:CC", "GO:MF", "REAC", "TF"),
  correction = "g_SCS",
  cutoff = 0.05,
  max_query_size = 10000,
  enrichment_cache = NULL
)
```

## Arguments

- basis_genes_list:

  A named list of character vectors, one gene set per factor (e.g. from
  splitting \[nmf_basis_genes()\]).

- background_genes:

  A character vector of background gene symbols.

- organism, sources, correction, cutoff:

  g:Profiler query settings; see \[NMFprofileR()\] for details.

- max_query_size:

  Integer; gene sets larger than this are skipped rather than sent to
  g:Profiler (default 10000).

- enrichment_cache:

  Optional path to a directory used to cache g:Profiler results by query
  hash; \`NULL\` (the default) disables caching. See \[NMFprofileR()\].

## Value

A list with two elements: \`per_factor\` (a list of \`gost\` result
objects, one per factor, \`NULL\` where there were no results) and
\`combined\` (a data frame of enrichment across the union of all factor
genes).

## See also

\[nmf_basis_genes()\]

## Examples

``` r
if (FALSE) { # \dontrun{
bg <- rownames(m)
gene_sets <- split(basis_genes$Gene, basis_genes$Factor)
enr <- nmf_enrichment(gene_sets, background_genes = bg)
} # }
```
