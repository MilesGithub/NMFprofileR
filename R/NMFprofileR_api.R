#' Preprocess an expression matrix for NMF
#'
#' Validates the input (numeric matrix, gene-symbol rownames, duplicate-symbol
#' resolution, all-zero row/column checks) and applies the standard filtering
#' (optional gene-list filter, mean-expression threshold, variance quantile,
#' non-negativity). This is the first stage of the [NMFprofileR()] pipeline,
#' exposed for composable use.
#'
#' @inheritParams NMFprofileR
#'
#' @return A filtered, non-negative numeric matrix ready for [nmf_fit()].
#' @seealso [nmf_fit()], [NMFprofileR()]
#' @export
#'
#' @examples
#' data("example_expression_data")
#' m <- nmf_preprocess(example_expression_data, expression_threshold = 0,
#'                     variance_quantile = 0)
#' dim(m)
nmf_preprocess <- function(expression_data,
                           expression_threshold = 10.0,
                           variance_quantile = 0.7,
                           gene_list_filter_file = NULL,
                           on_duplicate_genes = c("collapse_max", "collapse_mean", "error"),
                           variance_scale = c("raw", "log")) {
  on_duplicate_genes <- match.arg(on_duplicate_genes)
  variance_scale <- match.arg(variance_scale)
  expression_data <- validate_expression_input(expression_data, on_duplicate_genes = on_duplicate_genes)
  expr_matrix <- preprocess_matrix(expression_data, expression_threshold, variance_quantile,
                                   gene_list_filter_file, variance_scale = variance_scale)
  if (nrow(expr_matrix) < 2) {
    stop(
      "Fewer than 2 genes remain after filtering, so NMF cannot proceed. ",
      "Check `expression_threshold` (default 10 assumes non-log-scale input; ",
      "on log-scaled data it removes almost everything) and `variance_quantile`.",
      call. = FALSE
    )
  }
  expr_matrix
}

#' Fit consensus NMF at a single rank
#'
#' A thin, in-memory wrapper around [NMF::nmf()] with the package's platform-aware
#' parallel handling and error trapping. Returns the fit object without writing
#' anything to disk.
#'
#' @param expr_matrix A preprocessed, non-negative numeric matrix (see
#'   [nmf_preprocess()]).
#' @param rank The factorization rank (integer).
#' @param method The NMF algorithm (e.g. `"brunet"`).
#' @param nrun The number of consensus runs.
#' @param seed An integer seed for reproducibility.
#'
#' @return An `NMFfit` object, or `NULL` if the fit fails.
#' @seealso [nmf_preprocess()], [nmf_basis_genes()], [nmf_sample_assignments()]
#' @export
#'
#' @examples
#' \dontrun{
#' m <- nmf_preprocess(example_expression_data, expression_threshold = 0,
#'                     variance_quantile = 0)
#' fit <- nmf_fit(m, rank = 3, nrun = 10)
#' }
nmf_fit <- function(expr_matrix, rank, method = "brunet", nrun = 20, seed = 123456) {
  is_windows <- tolower(.Platform$OS.type) == "windows"
  nmf_options <- if (is_windows) {
    list(parallel = 0, verbose = TRUE)
  } else {
    list(parallel = min(nrun, max(1, tryCatch(parallel::detectCores(logical = TRUE), error = function(e) 1))),
         verbose = TRUE)
  }
  pbackend_arg <- if (is_windows) NA else NULL

  tryCatch({
    args <- list(x = expr_matrix, rank = rank, method = method, nrun = nrun, seed = seed, .options = nmf_options)
    if (!is.null(pbackend_arg)) args$.pbackend <- pbackend_arg
    do.call(NMF::nmf, args)
  }, error = function(e) {
    warning(sprintf("NMF fit failed for rank %s: %s", rank, conditionMessage(e)), call. = FALSE)
    NULL
  })
}

#' Assign genes to their dominant NMF factor
#'
#' Assigns every gene to its dominant basis factor (argmax over the basis matrix,
#' via [NMF::predict()]), returning a tidy table ordered by factor then gene.
#'
#' @param fit An `NMFfit` object from [nmf_fit()].
#'
#' @return A tibble with columns `Gene` and `Factor` (a factor with levels
#'   `Factor_1 ... Factor_k`).
#' @seealso [nmf_fit()], [nmf_enrichment()]
#' @export
#'
#' @examples
#' \dontrun{
#' basis_genes <- nmf_basis_genes(fit)
#' split(basis_genes$Gene, basis_genes$Factor)
#' }
nmf_basis_genes <- function(fit) {
  W <- NMF::basis(fit)
  k <- ncol(W)

  raw_gene_pred <- tryCatch(
    NMF::predict(fit, what = "features", dmatrix = TRUE),
    error = function(e) NULL
  )
  gene_preds <- normalize_predict(raw_gene_pred, nrow(W), axis = "features")

  df <- tibble::tibble(
    Gene = rownames(W),
    Factor = paste0("Factor_", gene_preds)
  )
  df$Factor <- factor(df$Factor, levels = paste0("Factor_", 1:k))
  df[order(df$Factor, df$Gene), ]
}

#' Extract factor-specific marker genes by specificity score
#'
#' Unlike [nmf_basis_genes()] (which assigns *every* gene to its dominant factor
#' by argmax), this selects only the genes that are specifically associated with
#' each factor, using `NMF::extractFeatures()` (the Kim & Park specificity
#' scoring). The result is a smaller, sharper set of marker genes per factor;
#' many genes are intentionally left unassigned.
#'
#' @param fit An `NMFfit` object from [nmf_fit()].
#'
#' @return A tibble with columns `Gene` and `Factor`. Factors with no specific
#'   markers contribute no rows.
#' @seealso [nmf_basis_genes()], [nmf_enrichment()]
#' @export
#'
#' @examples
#' \dontrun{
#' markers <- nmf_marker_genes(fit)
#' split(markers$Gene, markers$Factor)
#' }
nmf_marker_genes <- function(fit) {
  W <- NMF::basis(fit)
  gene_names <- rownames(W)
  k <- ncol(W)

  feats <- tryCatch(NMF::extractFeatures(fit), error = function(e) NULL)

  per_factor <- lapply(seq_len(k), function(i) {
    idx <- if (is.list(feats) && length(feats) >= i) feats[[i]] else NULL
    if (is.null(idx) || length(idx) == 0 || all(is.na(idx))) return(NULL)
    idx <- idx[!is.na(idx)]
    genes <- if (is.numeric(idx)) gene_names[idx] else as.character(idx)
    genes <- genes[!is.na(genes)]
    if (length(genes) == 0) return(NULL)
    tibble::tibble(Gene = genes, Factor = paste0("Factor_", i))
  })

  df <- do.call(rbind, per_factor)
  if (is.null(df)) df <- tibble::tibble(Gene = character(), Factor = character())
  df$Factor <- factor(df$Factor, levels = paste0("Factor_", 1:k))
  df[order(df$Factor, df$Gene), ]
}

#' Assign samples to their dominant NMF factor with silhouette widths
#'
#' Assigns every sample to its dominant coefficient factor (argmax over the
#' coefficient matrix) and attaches the consensus-matrix silhouette width for
#' each sample.
#'
#' @param fit An `NMFfit` object from [nmf_fit()].
#' @param verbose Logical; passed to the silhouette computation for diagnostic
#'   messages.
#'
#' @return A data frame with columns `SampleID`, `Dominant_Factor`, and
#'   `Silhouette_NMF`, ordered by factor then descending silhouette.
#' @seealso [nmf_fit()]
#' @export
#'
#' @examples
#' \dontrun{
#' nmf_sample_assignments(fit)
#' }
nmf_sample_assignments <- function(fit, verbose = FALSE) {
  H <- NMF::coef(fit)
  k <- nrow(H)

  sa <- tibble::tibble(
    SampleID = colnames(H),
    Dominant_Factor = paste0("Factor_", apply(H, 2, which.max))
  )
  sil <- compute_sample_silhouette(fit, H, k, verbose = verbose)
  sa <- merge(sa, sil, by = "SampleID", all.x = TRUE, sort = FALSE)
  sa$Dominant_Factor <- factor(sa$Dominant_Factor, levels = paste0("Factor_", 1:k))
  sa[order(sa$Dominant_Factor, -sa$Silhouette_NMF, sa$SampleID), ]
}

#' Run g:Profiler enrichment on NMF factor gene sets
#'
#' Runs per-factor and combined g:Profiler functional enrichment on a list of
#' factor gene sets, in memory (nothing is written to disk).
#'
#' @param basis_genes_list A named list of character vectors, one gene set per
#'   factor (e.g. from splitting [nmf_basis_genes()]).
#' @param background_genes A character vector of background gene symbols.
#' @param organism,sources,correction,cutoff g:Profiler query settings; see
#'   [NMFprofileR()] for details.
#'
#' @return A list with two elements: `per_factor` (a list of `gost` result
#'   objects, one per factor, `NULL` where there were no results) and `combined`
#'   (a data frame of enrichment across the union of all factor genes).
#' @seealso [nmf_basis_genes()]
#' @export
#'
#' @examples
#' \dontrun{
#' bg <- rownames(m)
#' gene_sets <- split(basis_genes$Gene, basis_genes$Factor)
#' enr <- nmf_enrichment(gene_sets, background_genes = bg)
#' }
nmf_enrichment <- function(basis_genes_list,
                           background_genes,
                           organism = "hsapiens",
                           sources = c("GO:BP", "GO:CC", "GO:MF", "REAC", "TF"),
                           correction = "g_SCS",
                           cutoff = 0.05) {
  per_factor <- perform_enrichment(
    basis_genes_list = basis_genes_list, organism = organism, cutoff = cutoff,
    correction = correction, background_genes = background_genes, sources = sources,
    output_dir = NULL, k = NA, write_files = FALSE
  )
  combined <- perform_combined_enrichment(
    basis_genes_list = basis_genes_list, organism = organism, cutoff = cutoff,
    correction = correction, background_genes = background_genes, sources = sources,
    output_dir = NULL, k = NA, write_files = FALSE
  )
  list(per_factor = per_factor, combined = combined)
}

#' Tidy NMF rank-selection diagnostics
#'
#' Lays the key rank-selection metrics side by side across the tested ranks so
#' the "authoritative" component count can be chosen (typically the largest
#' cophenetic correlation before it drops, or a silhouette peak). Works on the
#' `rank_metrics` produced by [NMFprofileR()] (via `NMF::nmfEstimateRank`).
#'
#' @param x An `nmf_profile` object from [NMFprofileR()], or a `rank_metrics`
#'   data frame directly.
#' @param metrics Character vector of metrics to extract; matched case-insensitively
#'   and by prefix (so `"silhouette"` matches `silhouette.consensus`). Metrics
#'   not present in the data are silently skipped.
#'
#' @return A tibble in long form with columns `rank`, `metric`, and `value`.
#' @export
#'
#' @examples
#' rm <- data.frame(rank = 2:4, cophenetic = c(0.99, 0.95, 0.88),
#'                  dispersion = c(0.9, 0.8, 0.7), silhouette = c(0.6, 0.5, 0.4))
#' nmf_rank_diagnostics(rm)
nmf_rank_diagnostics <- function(x, metrics = c("cophenetic", "dispersion", "silhouette")) {
  df <- if (inherits(x, "nmf_profile")) x$rank_metrics else as.data.frame(x)
  if (is.null(df) || nrow(df) == 0) {
    return(tibble::tibble(rank = numeric(), metric = character(), value = numeric()))
  }

  cols <- names(df)
  rank_hit <- cols[tolower(cols) == "rank"]
  ranks <- if (length(rank_hit) >= 1) df[[rank_hit[1]]] else seq_len(nrow(df))

  parts <- lapply(metrics, function(m) {
    hit <- cols[tolower(cols) == tolower(m)]
    if (length(hit) == 0) hit <- cols[startsWith(tolower(cols), tolower(m))]
    if (length(hit) == 0) return(NULL)
    tibble::tibble(rank = as.numeric(ranks), metric = m, value = as.numeric(df[[hit[1]]]))
  })
  dplyr::bind_rows(parts)
}

#' Print an NMFprofileR result
#'
#' @param x An `nmf_profile` object returned (invisibly) by [NMFprofileR()].
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.nmf_profile <- function(x, ...) {
  ranks <- names(x$fits)
  n_samples <- if (length(x$sample_assignments)) nrow(x$sample_assignments[[1]]) else NA_integer_
  ver <- x$provenance$gprofiler_version

  cat("<nmf_profile>\n")
  cat(sprintf("  Ranks fitted : %s\n", if (length(ranks)) paste(ranks, collapse = ", ") else "none"))
  cat(sprintf("  Samples      : %s\n", n_samples))
  cat(sprintf("  Runtime      : %s\n", format(x$runtime)))
  cat(sprintf("  g:Profiler   : %s\n", if (is.null(ver) || is.na(ver)) "not captured" else ver))
  if (is.null(x$output_dirs)) {
    cat("  Output       : in memory only (write_files = FALSE)\n")
  } else {
    cat(sprintf("  Output dir   : %s\n", x$output_dirs$main))
  }
  cat("  Access $fits, $basis_genes, $sample_assignments, $enrichment, $consolidated_summary_df\n")
  invisible(x)
}

#' Summarize an NMFprofileR result
#'
#' @param object An `nmf_profile` object from [NMFprofileR()].
#' @param ... Ignored.
#' @return The consolidated per-factor summary data frame across all ranks.
#' @export
summary.nmf_profile <- function(object, ...) {
  object$consolidated_summary_df
}