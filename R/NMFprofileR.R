#' Perform a Full Multi-Rank NMF and Functional Enrichment Analysis
#'
#' This function is the main entry point for the NMFprofileR workflow. It
#' orchestrates the entire analysis pipeline, from preprocessing the input
#' expression matrix to performing a multi-rank NMF analysis, extracting basis
#' genes, running functional enrichment, and generating a comprehensive suite of
#' visualizations and summary reports.
#'
#' @param expression_data A data frame or matrix of normalized, non-negative
#'   gene expression values. Rows should represent genes (with gene symbols as
#'   rownames) and columns should represent samples.
#' @param nmf_rank An integer vector specifying the factorization ranks (k) to
#'   test (e.g., `2:5`).
#' @param output_prefix A character string defining the base path and file name
#'   prefix for all output files and directories. For example, `"results/MyAnalysis"`.
#' @param gene_list_filter_file An optional character string path to a plain text
#'   file containing a list of gene symbols (one per line). If provided, the
#'   expression matrix will be pre-filtered to include only these genes.
#' @param expression_threshold A numeric value. Genes with a mean expression
#'   across all samples below this threshold will be removed prior to NMF.
#' @param variance_quantile A numeric value between 0 and 1. After mean-based
#'   filtering, genes with variance below this quantile are removed. This step
#'   enriches for more informative genes.
#' @param nmf_method A character string specifying the NMF algorithm to use.
#'   Defaults to 'brunet'. See the `NMF::nmf` documentation for other options.
#' @param nmf_nrun An integer specifying the number of runs for the consensus NMF.
#'   A higher number (e.g., 50-100) provides more stable results.
#' @param nmf_seed An integer to use as the random seed for reproducibility.
#' @param gprofiler_organism A character string specifying the organism for
#'   g:Profiler enrichment analysis (e.g., 'hsapiens', 'mmusculus').
#' @param gprofiler_sources A character vector of data sources to query in
#'   g:Profiler (e.g., c("GO:BP", "REAC", "KEGG")).
#' @param gprofiler_correction A character string specifying the multiple testing
#'   correction method used by g:Profiler.
#' @param gprofiler_cutoff A numeric value for the significance threshold for
#'   enrichment results.
#' @param enrichment_plot_top_n An integer specifying the number of top enriched
#'   terms to display in the summary dot plots.
#' @param verbose A logical value. If TRUE, prints detailed messages. Defaults to TRUE.
#'
#' @importFrom grDevices dev.control dev.off pdf recordPlot replayPlot
#'
#' @return A single data frame (`consolidated_summary_df`) containing a detailed
#'   summary for every factor from every rank tested. All detailed results, plots,
#'   and log files are saved to disk at the location specified by `output_prefix`.
#'
#' @export
#'
NMFprofileR <- function(
    expression_data,
    nmf_rank = 2:5,
    output_prefix,
    gene_list_filter_file = NULL,
    expression_threshold = 10.0,
    variance_quantile = 0.7,
    nmf_method = "brunet",
    nmf_nrun = 20,
    nmf_seed = 123456,
    gprofiler_organism = "hsapiens",
    gprofiler_sources = c("GO:BP", "GO:CC", "GO:MF", "REAC", "TF"),
    gprofiler_correction = "g_SCS",
    gprofiler_cutoff = 0.05,
    enrichment_plot_top_n = 50,
    verbose = FALSE
) {
  # --- Helper for verbose/debug messages ---
  vmsg <- function(...) {
    if (isTRUE(verbose)) cli::cli_alert_info(...)
  }

  time_start <- Sys.time()
  cli::cli_h1("Starting NMFprofileR Pipeline")

  if (!is.data.frame(expression_data) && !is.matrix(expression_data)) {
    stop("`expression_data` must be a data frame or matrix.", call. = FALSE)
  }

  # sanitize output_prefix (remove trailing slashes)
  output_prefix <- sub("/+$", "", as.character(output_prefix))
  if (nzchar(output_prefix) == FALSE) stop("`output_prefix` must be a non-empty string.", call. = FALSE)

  dirs <- setup_directories(output_prefix)
  file_prefix <- basename(output_prefix)

  # Open log sink safely; ensure all sinks closed on exit
  log_file_path <- file.path(dirs$main, paste0(file_prefix, "_run_log.txt"))
  sink(log_file_path, append = TRUE, split = TRUE)
  on.exit({
    # close any opened sinks
    while (sink.number() > 0) sink(NULL)
    cli::cli_alert_success("Log file saved to: {.path {log_file_path}}")
  }, add = TRUE)

  # --- Preprocessing ---
  cli::cli_h2("Step 1: Preprocessing Expression Matrix")
  vmsg("Converting and filtering expression matrix...")
  expr_matrix <- preprocess_matrix(expression_data, expression_threshold, variance_quantile, gene_list_filter_file)
  final_nmf_genes <- rownames(expr_matrix)
  cli::cli_alert_info("Final matrix for NMF: {nrow(expr_matrix)} genes x {ncol(expr_matrix)} samples.")

  # --- NMF Rank Estimation Survey (safe parallel handling) ---
  cli::cli_h2("Step 2: Performing NMF Rank Estimation Survey")
  is_windows <- tolower(.Platform$OS.type) == "windows"

  estim_real <- tryCatch({
    nmf_est_args <- list(
      x = expr_matrix,
      range = nmf_rank,
      method = nmf_method,
      nrun = nmf_nrun,
      seed = nmf_seed
    )

    if (is_windows) {
      nmf_est_args$.options <- list(parallel = 0)
      nmf_est_args$.pbackend <- NA
    } else {
      available_cores <- tryCatch(parallel::detectCores(logical = TRUE), error = function(e) 1)
      nmf_est_args$.options <- list(parallel = min(nmf_nrun, max(1, available_cores)))
    }

    do.call(NMF::nmfEstimateRank, nmf_est_args)
  }, error = function(e) {
    warning("NMF rank estimation failed: ", conditionMessage(e), call. = FALSE)
    return(NULL)
  })

  all_rank_metrics_df <- data.frame()
  if (!is.null(estim_real)) {
    plot_path <- file.path(dirs$plots, "02_Rank_Survey_Plot.pdf")
    cli::cli_alert_info("Saving rank survey plot to: {.path {basename(plot_path)}}")
    grDevices::pdf(plot_path, width = 10, height = 7)
    suppressWarnings(print(plot(estim_real)))
    grDevices::dev.off()
    all_rank_metrics_df <- as.data.frame(summary(estim_real))
  } else {
    cli::cli_alert_warning("Rank estimation returned NULL; continuing to run specified ranks.")
  }

  # --- Fit NMF for each rank and downstream analysis ---
  cli::cli_h2("Step 3: Performing Deep Analysis for Each Rank")
  all_summaries_list <- list()
  nmf_rds_paths <- list()

  for (k in nmf_rank) {

    cli::cli_rule(left = paste("Processing Rank k =", k))
    vmsg("Running NMF for k = {k} ...")
    nmf_result <- run_nmf_for_rank(k, expr_matrix, nmf_method, nmf_nrun, nmf_seed, dirs$nmf_core, file_prefix)

    if (!is.null(nmf_result)) {

        # save path (in case run_nmf_for_rank didn't)
        nmf_rds <- file.path(dirs$nmf_core, paste0("NMF_Result_Object_Rank_k", k, ".rds"))
        if (file.exists(nmf_rds)) nmf_rds_paths[[as.character(k)]] <- nmf_rds

        # Extract W and H
        W <- NMF::basis(nmf_result)
        H <- NMF::coef(nmf_result)

        # --- Assign genes robustly (handle vector, list($predict), or membership matrix) ---
        normalize_predict <- function(raw, expected_length, axis = c("features", "samples")) {
          axis <- match.arg(axis)
          if (is.list(raw) && "predict" %in% names(raw)) {
            preds <- as.integer(raw$predict)
          } else if (is.atomic(raw) && length(raw) == expected_length) {
            preds <- as.integer(raw)
          } else if (is.matrix(raw) && ncol(raw) >= 1 && nrow(raw) == expected_length) {
            preds <- apply(raw, 1, which.max)
          } else if (is.matrix(raw) && ncol(raw) >= 1 && ncol(raw) == expected_length) {
            # fallback if transposed membership matrix
            preds <- apply(raw, 2, which.max)
          } else {
            preds <- rep(NA_integer_, expected_length)
          }
          if (length(preds) != expected_length) {
            preds <- preds[seq_len(min(length(preds), expected_length))]
            if (length(preds) < expected_length) preds <- c(preds, rep(NA_integer_, expected_length - length(preds)))
          }
          return(as.integer(preds))
        }

        # gene predictions
        raw_gene_pred <- tryCatch(
          NMF::predict(nmf_result, what = "features", dmatrix = TRUE),
          error = function(e) {
            return(NULL)
          }
        )

        gene_preds <- normalize_predict(raw_gene_pred, nrow(W), axis = "features")

        basis_genes_df <- tibble::tibble(
          Gene = rownames(W),
          Factor = paste0("Factor_", gene_preds)
        )

        # Order by Factor first, then alphabetically within factor
        basis_genes_df$Factor <- factor(
          basis_genes_df$Factor,
          levels = paste0("Factor_", 1:k)
        )
        basis_genes_df <- basis_genes_df[order(basis_genes_df$Factor, basis_genes_df$Gene), ]

        # Ensure directory exists
        dir.create(dirs$basis_genes, showWarnings = FALSE, recursive = TRUE)
        readr::write_tsv(
          basis_genes_df,
          file.path(dirs$basis_genes, paste0("Basis_Genes_Rank_k", k, ".tsv"))
        )

        # Split into list
        basis_genes_list <- split(basis_genes_df$Gene, basis_genes_df$Factor)


        # --- Enrichment (use named args to avoid positional matching issues) ---
        gost_objects_list <- perform_enrichment(
          basis_genes_list = basis_genes_list,
          organism = gprofiler_organism,
          cutoff = gprofiler_cutoff,
          correction = gprofiler_correction,
          background_genes = final_nmf_genes,
          sources = gprofiler_sources,
          output_dir = dirs$enrichment,
          k = k
        )

        # filter NULLs before binding
        valid_gost_objects <- gost_objects_list[!sapply(gost_objects_list, is.null)]
        all_gprofiler_results_dfs <- lapply(valid_gost_objects, function(gost_obj) gost_obj$result)
        combined_gprofiler_df <- dplyr::bind_rows(all_gprofiler_results_dfs)

        combined_enrichment_results_df <- perform_combined_enrichment(
          basis_genes_list = basis_genes_list,
          organism = gprofiler_organism,
          cutoff = gprofiler_cutoff,
          correction = gprofiler_correction,
          background_genes = final_nmf_genes,
          sources = gprofiler_sources,
          output_dir = dirs$enrichment,
          k = k
        )
        if (is.null(combined_enrichment_results_df)) combined_enrichment_results_df <- data.frame()

        # --- Sample assignments ---
        raw_sample_pred <- tryCatch(
          NMF::predict(nmf_result, what = "samples", dmatrix = TRUE),
          error = function(e) {
            return(NULL)
          }
        )

        sample_preds <- normalize_predict(raw_sample_pred, ncol(H), axis = "samples")

        sample_assignments <- tibble::tibble(
          SampleID = colnames(H),
          Dominant_Factor = paste0("Factor_", sample_preds)
        )

        # Order by Factor first, then alphabetically within factor
        sample_assignments$Dominant_Factor <- factor(
          sample_assignments$Dominant_Factor,
          levels = paste0("Factor_", 1:k)
        )
        sample_assignments <- sample_assignments[order(sample_assignments$Dominant_Factor,
                                                       sample_assignments$SampleID), ]

        dir.create(dirs$sample_assignments, showWarnings = FALSE, recursive = TRUE)
        readr::write_tsv(
          sample_assignments,
          file.path(dirs$sample_assignments, paste0("Sample_Assignments_Rank_k", k, ".tsv"))
        )

        # --- Rank summary & plots ---
        k_summary_df <- generate_rank_summary(k, W, H, sample_assignments, basis_genes_list, combined_gprofiler_df)
        all_summaries_list[[as.character(k)]] <- k_summary_df

        k_plots_dir <- file.path(dirs$plots, paste0("Rank_k", k))
        dir.create(k_plots_dir, showWarnings = FALSE, recursive = TRUE)
        generate_rank_plots(
          k = k,
          nmf_result = nmf_result,
          sample_assignments = sample_assignments,
          combined_gprofiler_df = combined_gprofiler_df,
          combined_enrichment_results_df = combined_enrichment_results_df,
          expr_matrix = expr_matrix,
          basis_genes = basis_genes_list,
          k_plots_dir = k_plots_dir,
          file_prefix = file_prefix,
          nrun = nmf_nrun,
          top_n = enrichment_plot_top_n,
          gost_objects_list = gost_objects_list
        )
      }
    }
    # --- Finalize and Global Plots ---
    cli::cli_h2("Step 4: Finalizing and Generating Global Summary Plots")
    consolidated_summary_df <- if (length(all_summaries_list) > 0) dplyr::bind_rows(all_summaries_list) else tibble::tibble()
    dir.create(dirs$summaries, showWarnings = FALSE, recursive = TRUE)
    if (nrow(consolidated_summary_df) > 0) readr::write_tsv(consolidated_summary_df, file.path(dirs$summaries, "Consolidated_Summary.tsv"))

    # generate global plots (guard if summary empty)
    tryCatch({
      generate_global_plots(all_rank_metrics_df, consolidated_summary_df, nmf_rank, nmf_nrun, dirs, file_prefix)
    }, error = function(e) {
      cli::cli_alert_warning("generate_global_plots() failed: {conditionMessage(e)}")
    })

    time_end <- Sys.time()
    runtime <- difftime(time_end, time_start)
    cli::cli_alert_success("NMFprofileR Pipeline finished successfully in {format(runtime)}.")

    # Return a programmatically useful list while keeping consolidated_summary_df as main item
    result <- list(
      consolidated_summary_df = consolidated_summary_df,
      rank_metrics = all_rank_metrics_df,
      nmf_rds = nmf_rds_paths,
      output_dirs = dirs,
      runtime = runtime
    )

  invisible(result)
}
