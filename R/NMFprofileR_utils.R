#' Setup a Standardized Directory Structure for NMF Analysis Results
#'
#' This internal helper function takes a base output prefix and creates a
#' full, hierarchical directory structure for storing all results from the
#' NMFprofileR workflow. It is designed to be robust and will not produce
#' warnings if the directories already exist.
#'
#' @param output_prefix A character string defining the base path and file name
#'   prefix for the analysis. For example, "results/MyAnalysis".
#'
#' @return A named list where each element is the full path to one of the
#'   created subdirectories (e.g., `dirs$plots`, `dirs$summaries`).
#'
#' @noRd
setup_directories <- function(output_prefix) {
  output_prefix <- sub("/+$", "", as.character(output_prefix))
  file_prefix <- basename(output_prefix)
  base_dir <- dirname(output_prefix)
  if (base_dir == ".") base_dir <- getwd()

  main_results_dir <- file.path(base_dir, paste0(file_prefix, "_Results"))

  dirs <- list(
    main = main_results_dir,
    nmf_core = file.path(main_results_dir, "NMF_Core_Objects"),
    basis_genes = file.path(main_results_dir, "Basis_Genes"),
    enrichment = file.path(main_results_dir, "Enrichment_Results"),
    sample_assignments = file.path(main_results_dir, "Sample_Assignments"),
    plots = file.path(main_results_dir, "Plots"),
    summaries = file.path(main_results_dir, "Summaries")
  )

  invisible(lapply(dirs, dir.create, showWarnings = FALSE, recursive = TRUE))

  return(dirs)
}

#' Capture a plot as a reusable grob without rendering to a visible device
#'
#' This function evaluates a plotting expression, captures it using
#' recordPlot(), and returns the captured plot object for later use.
#' Useful for grid-based graphics (e.g., ComplexHeatmap, ggplot2) or
#' when saving multiple plots programmatically.
#'
#' @param plot_expression Expression that generates a plot (e.g., ggplot object or base R plot call)
#' @return A recorded plot object (can be replayed with replayPlot)
capture_plot <- function(plot_expression) {
  pdf(NULL);
  dev.control(displaylist = "enable");
  plot_expression;
  p <- recordPlot();
  dev.off()

  return(p)
}

#' Preprocess a Gene Expression Matrix for NMF Analysis
#'
#' This internal helper function applies a series of filtering steps to prepare
#' a raw expression matrix for NMF. The steps include an optional pre-filter
#' with a user-provided gene list, filtering by mean expression, and filtering
#' by variance. It also ensures the final matrix is non-negative.
#'
#' @param expression_data A data frame or matrix of gene expression values.
#' @param threshold Numeric value for mean expression filtering.
#' @param v_quantile Numeric value (0-1) for variance filtering.
#' @param filter_file A character string path to an optional gene list file.
#'
#' @return A numeric matrix that has been filtered and is ready for NMF.
#'
#' @noRd
preprocess_matrix <- function(expression_data, threshold, v_quantile, filter_file) {
  expr_matrix <- as.matrix(expression_data)

  if (any(!is.finite(expr_matrix))) {
    cli::cli_alert_warning("Non-finite values (NA, NaN, Inf) detected. Imputing with 0.")
    expr_matrix[!is.finite(expr_matrix)] <- 0
  }

  if (!is.null(filter_file) && nzchar(filter_file) && file.exists(filter_file)) {
    cli::cli_alert_info("Applying gene filter list from: {.path {basename(filter_file)}}")
    gene_filter_list <- readr::read_lines(filter_file)
    genes_to_keep <- intersect(rownames(expr_matrix), gene_filter_list)
    expr_matrix <- expr_matrix[genes_to_keep, , drop = FALSE]
  }

  cli::cli_alert_info("Filtering genes by mean expression < {threshold}...")
  if (nrow(expr_matrix) == 0) return(expr_matrix)
  expr_matrix <- expr_matrix[rowMeans(expr_matrix, na.rm = TRUE) >= threshold, , drop = FALSE]

  if (nrow(expr_matrix) > 1) {
    cli::cli_alert_info("Filtering genes by variance < {v_quantile} quantile...")
    gene_vars <- apply(expr_matrix, 1, stats::var, na.rm = TRUE)
    pos_vars <- gene_vars[gene_vars > 0]
    if (length(pos_vars) > 0) {
      variance_cutoff <- stats::quantile(pos_vars, probs = v_quantile, na.rm = TRUE)
      expr_matrix <- expr_matrix[gene_vars >= variance_cutoff, , drop = FALSE]
    } else {
      cli::cli_alert_info("No positive variances found; skipping variance filtering.")
    }
  }

  if (any(expr_matrix < 0, na.rm = TRUE)) {
    cli::cli_alert_warning("Negative values detected. Setting to 0 for NMF.")
    expr_matrix[expr_matrix < 0] <- 0
  }

  return(expr_matrix)
}

#' Run Consensus NMF for a Single Rank
#'
#' A robust wrapper for the `NMF::nmf()` execution for a single rank `k`. It
#' includes error handling and saves the full NMF result object to disk.
#'
#' @param k The rank (integer) for which to run the NMF analysis.
#' @param expr_matrix The pre-processed, non-negative numeric matrix.
#' @param method The NMF algorithm to use (e.g., "brunet").
#' @param nrun The number of runs for the consensus NMF.
#' @param seed An integer seed for reproducibility.
#' @param nmf_core_dir The path to the output directory for saving the NMF object.
#' @param file_prefix The base prefix for the output file name.
#'
#' @return The `NMFfit` object, or `NULL` if the NMF execution fails.
#'
#' @noRd
run_nmf_for_rank <- function(k, expr_matrix, method, nrun, seed, nmf_core_dir, file_prefix) {
  cli::cli_alert_info("Running Consensus NMF (k={k}, nrun={nrun})...")

  is_windows <- tolower(.Platform$OS.type) == "windows"
  nmf_options <- if (is_windows) list(parallel = 0, verbose = TRUE) else list(parallel = min(nrun, max(1, tryCatch(parallel::detectCores(logical = TRUE), error = function(e) 1))), verbose = TRUE)
  pbackend_arg <- if (is_windows) NA else NULL

  nmf_result <- tryCatch({
    args <- list(x = expr_matrix, rank = k, method = method, nrun = nrun, seed = seed, .options = nmf_options)
    if (!is.null(pbackend_arg)) args$.pbackend <- pbackend_arg
    do.call(NMF::nmf, args)
  }, error = function(e) {
    warning(sprintf("NMF fit failed for rank %d: %s", k, conditionMessage(e)), call. = FALSE)
    return(NULL)
  })

  if (!is.null(nmf_result)) {
    save_path <- file.path(nmf_core_dir, paste0("NMF_Result_Object_Rank_k", k, ".rds"))
    dir.create(nmf_core_dir, showWarnings = FALSE, recursive = TRUE)
    saveRDS(nmf_result, save_path)
  }

  return(nmf_result)
}

#' Perform Functional Enrichment for a List of Gene Sets
#'
#' This function iterates over a list of gene sets (one for each NMF factor)
#' and performs g:Profiler functional enrichment analysis.
#'
#' @param basis_genes_list A named list where each element is a character vector of gene symbols.
#' @param organism A character string for the g:Profiler organism.
#' @param cutoff A numeric significance threshold.
#' @param correction A character string for the multiple testing correction method.
#' @param background_genes A character vector of background genes.
#' @param sources A character vector of g:Profiler data sources.
#' @param output_dir Path to save the enrichment results.
#' @param k The current rank, used for file naming.
#'
#' @return A list of data frames, where each data frame contains the enrichment
#'   results for a single factor. Returns `NULL` for factors with no results.
#'
#' @noRd
perform_enrichment <- function(basis_genes_list,
                               organism,
                               cutoff,
                               correction,
                               background_genes,
                               sources,
                               output_dir,
                               k) {
  cli::cli_alert_info("Performing per-factor functional enrichment (k={k})...")

  all_results <- lapply(names(basis_genes_list), function(factor_name) {
    gene_set <- basis_genes_list[[factor_name]]
    if (length(gene_set) < 5) return(NULL)

    gost_res <- tryCatch({

      named_query <- list(gene_set)
      names(named_query) <- factor_name
      suppressMessages(
        suppressWarnings(
          gprofiler2::gost(
            query = named_query,
            organism = organism,
            sources = sources,
            correction_method = correction,
            user_threshold = cutoff,
            custom_bg = background_genes,
            significant = TRUE
          )))
    }, error = function(e) {
      warning(sprintf("g:Profiler query failed for %s (k=%s): %s", factor_name, k, e$message), call. = FALSE)
      return(NULL)
    })

    if (!is.null(gost_res) && is.data.frame(gost_res$result) && nrow(gost_res$result) > 0) {
      readr::write_tsv(gost_res$result, file.path(output_dir, paste0("Enrichment_Rank_k", k, "_", factor_name, ".tsv")))
      return(gost_res)
    } else {
      return(NULL)
    }
  })

  return(all_results)
}


#' Perform Functional Enrichment on All Basis Genes Combined
#'
#' Runs a single g:Profiler query on the unique set of all basis genes from all factors.
#'
#' @param basis_genes_list A named list of gene sets.
#' @param organism A character string for the g:Profiler organism.
#' @param cutoff A numeric significance threshold.
#' @param correction A character string for the multiple testing correction method.
#' @param background_genes A character vector of background genes.
#' @param sources A character vector of g:Profiler data sources.
#' @param output_dir Path to save the enrichment results.
#' @param k The current rank, used for file naming.
#'
#' @return A data frame of enrichment results, or an empty data frame on failure.
#'
#' @noRd
perform_combined_enrichment <- function(basis_genes_list,
                                        organism,
                                        cutoff,
                                        correction,
                                        background_genes,
                                        sources,
                                        output_dir,
                                        k) {
  cli::cli_alert_info("Performing combined enrichment on all basis genes (k={k})...")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  combined_genes <- unique(unlist(basis_genes_list))
  if (length(combined_genes) < 5) return(data.frame())

  gost_res <- tryCatch({
    suppressMessages(
      suppressWarnings(
        gprofiler2::gost(
          query = combined_genes,
          organism = organism,
          sources = sources,
          correction_method = correction,
          user_threshold = cutoff,
          custom_bg = background_genes,
          significant = TRUE
    )))
  }, error = function(e) {
    warning(sprintf("Combined g:Profiler query failed for rank %d: %s", k, conditionMessage(e)), call. = FALSE)
    return(NULL)
  })

  if (!is.null(gost_res) && is.data.frame(gost_res$result) && nrow(gost_res$result) > 0) {
    out_path <- file.path(output_dir, paste0("Enrichment_Rank_k", k, "_All_Factors_Combined.tsv"))
    readr::write_tsv(gost_res$result, out_path)
    return(gost_res$result)
  }

  return(data.frame())
}

#' Generate a Summary Data Frame for a Single NMF Rank
#'
#' This internal function calculates key statistics for each factor within a single NMF
#' rank, including sample counts, gene counts, sparsity, and enrichment counts.
#' The enrichment source columns in the output will be ordered according to the
#' `gprofiler_sources` parameter.
#'
#' @param k The current rank (integer).
#' @param W The basis matrix (features x rank).
#' @param H The coefficient matrix (rank x samples).
#' @param sample_assignments A data frame mapping samples to factors.
#' @param basis_genes_list A list of basis genes for each factor.
#' @param combined_gprofiler_df A data frame of all per-factor enrichment results
#'   from gprofiler2.
#' @param gprofiler_sources A character vector of the g:Profiler sources that
#'   were queried, used to order the final output columns.
#'
#' @return A data frame summarizing the results for rank `k`.
#'
#' @importFrom tidyselect where all_of
#' @noRd
generate_rank_summary <- function(k, W, H, sample_assignments, basis_genes_list,
                                  combined_gprofiler_df, gprofiler_sources) {

  # Process each factor individually and store results in a list
  summary_list <- lapply(1:k, function(i) {
    f_name <- paste0("Factor_", i)

    # Count the number of significant terms per source for the current factor.
    # This correctly filters by the 'query' column from the gprofiler2 result.
    enrich_counts <- if (nrow(combined_gprofiler_df) > 0 && "query" %in% names(combined_gprofiler_df)) {
      combined_gprofiler_df %>%
        dplyr::filter(.data$query == f_name) %>%
        dplyr::count(.data$source) %>%
        tidyr::pivot_wider(names_from = "source", values_from = "n", values_fill = 0)
    } else {
      tibble::tibble() # Return an empty tibble if no enrichment data exists
    }

    # Assemble the main summary statistics for the factor
    summary_row <- tibble::tibble(
      Factor = f_name,
      Num_Samples = sum(sample_assignments$Dominant_Factor == f_name, na.rm = TRUE),
      Percent_Samples = 100 * .data$Num_Samples / max(1, nrow(sample_assignments)),
      Num_Basis_Genes = length(basis_genes_list[[f_name]]),
      Basis_Sparsity = NMF::sparseness(W[, i]),
      Coef_Sparsity = NMF::sparseness(H[i, ])
    )

    # Combine the main stats with the enrichment counts, if any were found
    if (nrow(enrich_counts) > 0) {
      dplyr::bind_cols(summary_row, enrich_counts)
    } else {
      summary_row
    }
  })

  # Combine summaries for all factors into a single data frame
  combined_df <- dplyr::bind_rows(summary_list)

  # Identify which of the user-provided source columns actually exist in our results.
  # This handles cases where some sources might not have had any significant terms.
  existing_source_cols <- intersect(gprofiler_sources, names(combined_df))

  # Define the full, ordered set of columns based on user input
  final_col_order <- c(
    "Factor", "Num_Samples", "Percent_Samples", "Num_Basis_Genes",
    "Basis_Sparsity", "Coef_Sparsity",
    existing_source_cols
  )

  # Select and reorder the columns. `all_of()` ensures this works even if
  # some expected columns are missing (e.g., no enrichment results at all).
  # We also add any other unexpected columns to the end to be safe.
  final_df <- combined_df %>%
    dplyr::select(tidyselect::all_of(final_col_order), dplyr::everything())

  # Add the Rank column at the beginning and clean up NAs in numeric columns
  final_df %>%
    dplyr::mutate(Rank = k, .before = 1) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~tidyr::replace_na(., 0)))
}


#' Compute per-sample silhouette width from NMF consensus matrix
#'
#' This function calculates sample-level silhouette widths using the consensus
#' matrix from an NMF run (with nrun > 1). It mirrors the silhouette values
#' shown in NMF consensus maps.
#'
#' @param nmf_result An NMFfit object (result of NMF::nmf).
#' @param H The NMF coefficient matrix (factors x samples).
#' @param k Integer, the rank of the factorization.
#' @param verbose Logical, whether to print diagnostic messages.
#'
#' @return A data frame with columns:
#' \itemize{
#'   \item SampleID – sample names (colnames of H).
#'   \item Silhouette_NMF – silhouette width for each sample (NA if unavailable).
#' }
#'
#' @importFrom NMF consensus predict
#' @importFrom stats as.dist dist
#' @importFrom cluster silhouette
#'
#' @keywords internal
compute_sample_silhouette <- function(nmf_result, H, k, verbose = FALSE) {
  tryCatch({
    # Consensus matrix
    cons_mat <- tryCatch(NMF::consensus(nmf_result), error = function(e) NULL)
    if (is.null(cons_mat) || !is.matrix(cons_mat) || nrow(cons_mat) < 2) {
      if (verbose) message(sprintf("[k=%d] Consensus matrix not available. Returning NA.", k))
      return(data.frame(SampleID = colnames(H), Silhouette_NMF = NA_real_))
    }

    # Assign row/col names if missing
    if (is.null(colnames(cons_mat)) || is.null(rownames(cons_mat))) {
      colnames(cons_mat) <- rownames(cons_mat) <- colnames(H)
    }

    # Distance from consensus
    dist_cons <- stats::as.dist(1 - cons_mat)

    # Predicted sample clusters
    sample_labels <- NMF::predict(nmf_result, what = "samples")
    if (is.factor(sample_labels) || is.character(sample_labels)) {
      parsed <- suppressWarnings(as.integer(gsub(".*?(\\d+)$", "\\1", as.character(sample_labels))))
      if (all(is.finite(parsed))) {
        sample_labels_int <- parsed
      } else {
        sample_labels_int <- as.integer(as.factor(sample_labels))
      }
    } else {
      sample_labels_int <- as.integer(sample_labels)
    }

    # Align length
    if (length(sample_labels_int) != attr(dist_cons, "Size")) {
      if (verbose) message(sprintf("[k=%d] Label length mismatch. Returning NA.", k))
      return(data.frame(SampleID = colnames(H), Silhouette_NMF = NA_real_))
    }

    # At least 2 clusters required
    if (length(unique(sample_labels_int)) < 2) {
      if (verbose) message(sprintf("[k=%d] Only one cluster present. Returning NA.", k))
      return(data.frame(SampleID = colnames(H), Silhouette_NMF = NA_real_))
    }

    # Compute silhouette
    sil_obj <- cluster::silhouette(sample_labels_int, dist_cons)
    sil_w <- as.numeric(sil_obj[, "sil_width"])
    sample_ids <- rownames(sil_obj)
    if (is.null(sample_ids) || length(sample_ids) == 0) sample_ids <- colnames(H)

    data.frame(SampleID = sample_ids, Silhouette_NMF = sil_w, stringsAsFactors = FALSE)
  }, error = function(e) {
    if (verbose) message(sprintf("[k=%d] Silhouette computation failed: %s", k, e$message))
    data.frame(SampleID = colnames(H), Silhouette_NMF = NA_real_)
  })
}
