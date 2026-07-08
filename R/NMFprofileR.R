#' Perform a Full Multi-Rank NMF and Functional Enrichment Analysis
#'
#' This function is the main entry point for the NMFprofileR workflow. It
#' orchestrates the entire analysis pipeline, from preprocessing the input
#' expression matrix to performing a multi-rank NMF analysis, extracting basis
#' genes, running functional enrichment, and generating a comprehensive suite of
#' visualizations and summary reports.
#'
#' @section Reproducibility:
#' For a fixed `nmf_seed` the NMF factorization is deterministic and yields
#' identical basis and coefficient matrices whether NMF runs sequentially or
#' across a parallel backend, within a fixed computational environment. Exact
#' numerical reproducibility across machines can still be affected by the BLAS
#' implementation and threading. Every run also writes a manifest and session
#' information (surfaced in the `provenance` element of the return value) so a
#' result can be traced to the parameters and g:Profiler database snapshot that
#' produced it.
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
#' @param verbose A logical value. If TRUE, prints detailed messages. Defaults to FALSE.
#' @param write_files A logical value. If TRUE (the default), all results, plots,
#'   and log files are written to disk under `output_prefix`. If FALSE, nothing
#'   is written and no directories are created; the pipeline computes in memory
#'   and returns its results only. Useful for programmatic use, examples, and
#'   tests (which must not write outside `tempdir()`).
#' @param on_duplicate_genes How to handle duplicated gene symbols in the
#'   rownames of `expression_data`: `"collapse_max"` (the default) keeps the
#'   per-sample maximum across duplicate rows (the most expressed probe),
#'   `"collapse_mean"` averages them, and `"error"` stops with a message.
#'
#' @importFrom grDevices dev.control dev.off pdf recordPlot replayPlot
#' @importFrom cluster silhouette
#'
#' @return Invisibly, an S3 object of class `nmf_profile` (a list underneath, so
#'   `$` access works) with `print()` and `summary()` methods. Elements:
#'   \describe{
#'     \item{`consolidated_summary_df`}{A data frame summarizing every factor
#'       from every rank tested.}
#'     \item{`rank_metrics`}{A data frame of quality metrics from the
#'       `NMF::nmfEstimateRank` rank survey (empty if the survey failed).}
#'     \item{`fits`}{A named list of the fitted `NMFfit` objects, keyed by rank.}
#'     \item{`basis_genes`}{A named list of per-rank basis-gene tables
#'       (gene to dominant factor).}
#'     \item{`sample_assignments`}{A named list of per-rank sample-assignment
#'       tables (sample to dominant factor, with silhouette widths).}
#'     \item{`enrichment`}{A list with `per_factor` and `combined` named lists of
#'       per-rank g:Profiler enrichment data frames.}
#'     \item{`nmf_rds`}{A named list of file paths to the saved NMF fit objects
#'       (empty when `write_files = FALSE`).}
#'     \item{`output_dirs`}{A named list of the output directories created for
#'       the run, or `NULL` when `write_files = FALSE`.}
#'     \item{`runtime`}{A `difftime` giving the total pipeline runtime.}
#'     \item{`provenance`}{A named list recording the captured g:Profiler
#'       database version (`gprofiler_version`), the `run_parameters` used, and
#'       the paths of the manifest and session-information files written
#'       (`files`), for reproducibility.}
#'   }
#'   Unless `write_files = FALSE`, all detailed results, plots, and log files are
#'   additionally written to disk under the location specified by `output_prefix`.
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
    verbose = FALSE,
    write_files = TRUE,
    on_duplicate_genes = c("collapse_max", "collapse_mean", "error")
) {
  on_duplicate_genes <- match.arg(on_duplicate_genes)
  # --- Helper for verbose/debug messages ---
  vmsg <- function(...) {
    if (isTRUE(verbose)) cli::cli_alert_info(...)
  }

  time_start <- Sys.time()
  cli::cli_h1("Starting NMFprofileR Pipeline")

  # --- Preprocessing (validates input, then filters; fails fast before any
  #     directories are created) ---
  cli::cli_h2("Step 1: Preprocessing Expression Matrix")
  vmsg("Validating and filtering expression matrix...")
  expr_matrix <- nmf_preprocess(
    expression_data,
    expression_threshold = expression_threshold,
    variance_quantile = variance_quantile,
    gene_list_filter_file = gene_list_filter_file,
    on_duplicate_genes = on_duplicate_genes
  )
  final_nmf_genes <- rownames(expr_matrix)
  cli::cli_alert_info("Final matrix for NMF: {nrow(expr_matrix)} genes x {ncol(expr_matrix)} samples.")

  # sanitize output_prefix (remove trailing slashes)
  output_prefix <- sub("/+$", "", as.character(output_prefix))
  if (nzchar(output_prefix) == FALSE) stop("`output_prefix` must be a non-empty string.", call. = FALSE)

  dirs <- setup_directories(output_prefix, create = write_files)
  file_prefix <- basename(output_prefix)

  # File-connection logger. This replaces a global output sink so the run no
  # longer hijacks the console (and no longer captures unrelated packages'
  # output). Console output continues via cli; milestone lines are appended to
  # the run log only when file writing is enabled.
  log_file_path <- file.path(dirs$main, paste0(file_prefix, "_run_log.txt"))
  log_con <- NULL
  if (isTRUE(write_files)) {
    dir.create(dirs$main, showWarnings = FALSE, recursive = TRUE)
    log_con <- file(log_file_path, open = "at")
    on.exit({
      if (!is.null(log_con)) try(close(log_con), silent = TRUE)
      cli::cli_alert_success("Log file saved to: {.path {log_file_path}}")
    }, add = TRUE)
  }
  log_line <- function(...) {
    if (is.null(log_con)) return(invisible())
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", paste0(...), "\n",
        sep = "", file = log_con)
  }
  log_line("NMFprofileR run started (seed=", nmf_seed, ", method=", nmf_method,
           ", ranks=", paste(range(nmf_rank), collapse = "-"), ").")
  log_line("Preprocessed matrix: ", nrow(expr_matrix), " genes x ", ncol(expr_matrix), " samples.")

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
    all_rank_metrics_df <- as.data.frame(summary(estim_real))
    if (isTRUE(write_files)) {
      plot_path <- file.path(dirs$plots, "02_Rank_Survey_Plot.pdf")
      cli::cli_alert_info("Saving rank survey plot to: {.path {basename(plot_path)}}")
      grDevices::pdf(plot_path, width = 10, height = 7)
      suppressWarnings(print(plot(estim_real)))
      grDevices::dev.off()
    }
  } else {
    cli::cli_alert_warning("Rank estimation returned NULL; continuing to run specified ranks.")
  }

  # --- Fit NMF for each rank and downstream analysis ---
  cli::cli_h2("Step 3: Performing Deep Analysis for Each Rank")
  all_summaries_list <- list()
  nmf_rds_paths <- list()
  gprofiler_version <- NA_character_  # captured from the first enrichment result

  # In-memory accumulators (returned so callers need not re-read files)
  fits_by_rank <- list()
  basis_genes_by_rank <- list()
  sample_assignments_by_rank <- list()
  enrichment_per_factor_by_rank <- list()
  enrichment_combined_by_rank <- list()

  for (k in nmf_rank) {

    cli::cli_rule(left = paste("Processing Rank k =", k))
    vmsg("Running NMF for k = {k} ...")
    nmf_result <- run_nmf_for_rank(k, expr_matrix, nmf_method, nmf_nrun, nmf_seed, dirs$nmf_core, file_prefix,
                                   write_files = write_files)

    if (!is.null(nmf_result)) {

        fits_by_rank[[as.character(k)]] <- nmf_result

        # save path (in case run_nmf_for_rank didn't)
        nmf_rds <- file.path(dirs$nmf_core, paste0("NMF_Result_Object_Rank_k", k, ".rds"))
        if (file.exists(nmf_rds)) nmf_rds_paths[[as.character(k)]] <- nmf_rds

        # Extract W and H
        W <- NMF::basis(nmf_result)
        H <- NMF::coef(nmf_result)

        # --- Assign genes to their dominant factor (composable stage) ---
        basis_genes_df <- nmf_basis_genes(nmf_result)
        basis_genes_by_rank[[as.character(k)]] <- basis_genes_df

        if (isTRUE(write_files)) {
          dir.create(dirs$basis_genes, showWarnings = FALSE, recursive = TRUE)
          readr::write_tsv(
            basis_genes_df,
            file.path(dirs$basis_genes, paste0("Basis_Genes_Rank_k", k, ".tsv"))
          )
        }

        # Split into list
        basis_genes_list <- split(basis_genes_df$Gene, basis_genes_df$Factor)
        log_line("Rank k=", k, ": fitted; ", nrow(W), " genes assigned to ", k, " factors.")


        # --- Enrichment (use named args to avoid positional matching issues) ---
        gost_objects_list <- perform_enrichment(
          basis_genes_list = basis_genes_list,
          organism = gprofiler_organism,
          cutoff = gprofiler_cutoff,
          correction = gprofiler_correction,
          background_genes = final_nmf_genes,
          sources = gprofiler_sources,
          output_dir = dirs$enrichment,
          k = k,
          write_files = write_files
        )

        # filter NULLs before binding
        valid_gost_objects <- gost_objects_list[!sapply(gost_objects_list, is.null)]
        all_gprofiler_results_dfs <- lapply(valid_gost_objects, function(gost_obj) gost_obj$result)
        combined_gprofiler_df <- dplyr::bind_rows(all_gprofiler_results_dfs)

        # Capture the g:Profiler database version once, for provenance
        if (is.na(gprofiler_version)) {
          gprofiler_version <- extract_gprofiler_version(valid_gost_objects)
        }

        combined_enrichment_results_df <- perform_combined_enrichment(
          basis_genes_list = basis_genes_list,
          organism = gprofiler_organism,
          cutoff = gprofiler_cutoff,
          correction = gprofiler_correction,
          background_genes = final_nmf_genes,
          sources = gprofiler_sources,
          output_dir = dirs$enrichment,
          k = k,
          write_files = write_files
        )
        if (is.null(combined_enrichment_results_df)) combined_enrichment_results_df <- data.frame()
        enrichment_per_factor_by_rank[[as.character(k)]] <- combined_gprofiler_df
        enrichment_combined_by_rank[[as.character(k)]] <- combined_enrichment_results_df

        # --- Sample assignments + silhouette (composable stage) ---
        sample_assignments <- nmf_sample_assignments(nmf_result, verbose = TRUE)
        sample_assignments_by_rank[[as.character(k)]] <- sample_assignments

        # Save
        if (isTRUE(write_files)) {
          readr::write_tsv(
            sample_assignments,
            file.path(dirs$sample_assignments, paste0("Sample_Assignments_Rank_k", k, ".tsv"))
          )
        }

        # --- Rank summary & plots ---
        k_summary_df <- generate_rank_summary(
          k = k,
          W = W,
          H = H,
          sample_assignments = sample_assignments,
          basis_genes_list = basis_genes_list,
          combined_gprofiler_df = combined_gprofiler_df,
          gprofiler_sources = gprofiler_sources
        )
        all_summaries_list[[as.character(k)]] <- k_summary_df

        if (isTRUE(write_files)) {
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
            gost_objects_list = gost_objects_list,
            nmf_seed = nmf_seed
          )
        }
      }
    }
    # --- Finalize and Global Plots ---
    cli::cli_h2("Step 4: Finalizing and Generating Global Summary Plots")
    consolidated_summary_df <- if (length(all_summaries_list) > 0) dplyr::bind_rows(all_summaries_list) else tibble::tibble()
    if (isTRUE(write_files)) {
      dir.create(dirs$summaries, showWarnings = FALSE, recursive = TRUE)
      if (nrow(consolidated_summary_df) > 0) readr::write_tsv(consolidated_summary_df, file.path(dirs$summaries, "Consolidated_Summary.tsv"))

      # generate global plots (guard if summary empty)
      tryCatch({
        generate_global_plots(all_rank_metrics_df, consolidated_summary_df, nmf_rank, nmf_nrun, dirs, file_prefix)
      }, error = function(e) {
        cli::cli_alert_warning("generate_global_plots() failed: {conditionMessage(e)}")
      })
    }

    time_end <- Sys.time()
    runtime <- difftime(time_end, time_start)
    cli::cli_alert_success("NMFprofileR Pipeline finished successfully in {format(runtime)}.")

    # --- Provenance: record the exact parameters, g:Profiler snapshot, and
    #     session information that produced this run (reproducibility hardening).
    run_parameters <- list(
      nmf_rank             = nmf_rank,
      nmf_method           = nmf_method,
      nmf_nrun             = nmf_nrun,
      nmf_seed             = nmf_seed,
      expression_threshold = expression_threshold,
      variance_quantile    = variance_quantile,
      gprofiler_organism   = gprofiler_organism,
      gprofiler_sources    = gprofiler_sources,
      gprofiler_correction = gprofiler_correction,
      gprofiler_cutoff     = gprofiler_cutoff
    )
    provenance <- list(
      gprofiler_version = gprofiler_version,
      run_parameters    = run_parameters,
      files             = NULL
    )
    if (isTRUE(write_files)) {
      provenance$files <- tryCatch(
        write_run_provenance(dirs, file_prefix, run_parameters, gprofiler_version, runtime),
        error = function(e) {
          cli::cli_alert_warning("Failed to write run provenance: {conditionMessage(e)}")
          NULL
        }
      )
    }
    log_line("Run finished in ", format(runtime), "; ",
             length(all_summaries_list), " rank(s) completed.")

    # Assemble the in-memory result. Keeps all previously documented elements
    # (so `$consolidated_summary_df` etc. still work) and adds the fitted
    # objects, assignments, and enrichment tables so callers need not re-read
    # files. Returned as an S3 `nmf_profile` (a list underneath).
    result <- list(
      consolidated_summary_df = consolidated_summary_df,
      rank_metrics = all_rank_metrics_df,
      fits = fits_by_rank,
      basis_genes = basis_genes_by_rank,
      sample_assignments = sample_assignments_by_rank,
      enrichment = list(
        per_factor = enrichment_per_factor_by_rank,
        combined   = enrichment_combined_by_rank
      ),
      nmf_rds = nmf_rds_paths,
      output_dirs = if (isTRUE(write_files)) dirs else NULL,
      runtime = runtime,
      provenance = provenance
    )
    class(result) <- "nmf_profile"

  invisible(result)
}
