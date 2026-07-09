#' Perform a Full Multi-Rank NMF and Functional Enrichment Analysis
#'
#' This function is the main entry point for the NMFprofileR workflow. It
#' orchestrates the entire analysis pipeline, from preprocessing the input
#' expression matrix to performing a multi-rank NMF analysis, extracting basis
#' genes, running functional enrichment, and generating a comprehensive suite of
#' visualizations and summary reports.
#'
#' @section Reproducibility:
#' NMF is run sequentially (not across a parallel backend), which is both what
#' keeps results reproducible for a fixed `nmf_seed` and what avoids a failure
#' mode of NMF's parallel execution when called from within a package. For a
#' fixed `nmf_seed` the factorization is therefore deterministic within a fixed
#' computational environment; exact numerical reproducibility across machines
#' can still be affected by the BLAS implementation and threading. Every run
#' also writes a manifest and session information (surfaced in the `provenance`
#' element of the return value) so a result can be traced to the parameters and
#' g:Profiler database snapshot that produced it.
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
#' @param variance_scale Either `"raw"` (default) or `"log"`. The scale on which
#'   gene variance is measured for the highly-variable-gene filter. `expression_data`
#'   is expected to be normalized, non-negative values (e.g. normalized counts or
#'   TPM); on such data `expression_threshold = 10` and raw-scale variance are
#'   meaningful. `"log"` ranks variability on `log2(x + 1)`, reducing the
#'   mean-variance confound of count data. The matrix fed to NMF is always the
#'   original (raw, non-negative) scale regardless of this setting.
#' @param nmf_method A character string specifying the NMF algorithm to use.
#'   Defaults to 'brunet'. See the `NMF::nmf` documentation for other options.
#' @param nmf_nrun An integer specifying the number of runs for the consensus NMF.
#'   A higher number (e.g., 50-100) provides more stable results.
#' @param nmf_seed An integer to use as the random seed for reproducibility.
#' @param gprofiler_organism A character string specifying the organism for
#'   g:Profiler enrichment analysis (e.g., 'hsapiens', 'mmusculus').
#' @param gprofiler_sources A character vector of data sources to query in
#'   g:Profiler (e.g., c("GO:BP", "REAC", "KEGG")).
#' @param gprofiler_correction A character string specifying the multiple-testing
#'   correction method used by g:Profiler. Defaults to `"g_SCS"`, g:Profiler's
#'   native Set Counts and Sizes method, which accounts for the overlapping,
#'   hierarchical structure of GO and pathway terms and is generally more
#'   appropriate (and less conservative) there than Bonferroni or Benjamini-
#'   Hochberg FDR. The query uses a custom background of the genes that survive
#'   preprocessing (`custom_bg`), so enrichment is assessed against the tested
#'   gene universe rather than the whole genome. Alternatives (`"fdr"`,
#'   `"bonferroni"`) are accepted for reviewer familiarity at some loss of power.
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
#' @param emit_marker_genes A logical value. If TRUE, in addition to the argmax
#'   basis-gene assignment the pipeline also extracts factor-specific marker
#'   genes (`NMF::extractFeatures`, Kim-Park specificity) and runs enrichment on
#'   them, emitting `Marker_Genes_*` and `Marker_Enrichment_*` outputs alongside
#'   the argmax results. Defaults to FALSE: markers are opt-in because computing
#'   them runs a second g:Profiler query per rank (roughly doubling the network
#'   work), which matters most for batch use. The argmax assignment is always the
#'   primary output. Regardless of this setting the return value always carries
#'   the `marker_genes` and `enrichment$markers` elements (empty when FALSE).
#' @param custom_theme An optional `ggplot2` theme object applied to the
#'   ggplot-based plots (the factor-summary bar plot, enrichment dot plots, and
#'   the sample UMAP) in place of the package's built-in theme. `NULL` (the
#'   default) uses the built-in theme.
#' @param factor_palette An optional character vector of colours used to colour
#'   NMF factors across the plots. `NULL` (the default) uses the built-in
#'   palette. When a rank has more factors than supplied colours, the palette is
#'   extended with `grDevices::hcl.colors()`.
#' @param umap_n_neighbors Integer; the number of neighbours for the
#'   sample-coefficient UMAP embedding (default 15). The UMAP is only drawn when
#'   the cohort has more samples than this value, and the effective neighbour
#'   count is capped at `n_samples - 1` so smaller cohorts still embed.
#' @param run_id An optional identifier for this run. When supplied it is stamped
#'   as a leading `Run_ID` column in `consolidated_summary_df` (and the on-disk
#'   `Consolidated_Summary.tsv`) and recorded in the run manifest, so results
#'   from many runs can be pooled and traced back to their source.
#' @param nmf_parallel A logical value. If `FALSE` (the default) the per-rank
#'   consensus runs are fitted sequentially. If `TRUE` they are spread across a
#'   local cluster; see [nmf_fit()] for details. For a fixed `nmf_seed` the
#'   parallel path yields the same factorization as the sequential one. The
#'   rank-estimation survey always runs sequentially.
#' @param nmf_cores Number of worker processes to use when
#'   `nmf_parallel = TRUE`. `NULL` (the default) uses one fewer than the number
#'   of detected cores.
#' @param max_query_size An integer. Factor gene sets larger than this are
#'   skipped (with a warning) rather than sent to g:Profiler, guarding against
#'   pathological oversized queries. Defaults to 10000.
#' @param enrichment_cache An optional path to a directory used to cache
#'   g:Profiler results. When set, each query is keyed by a hash of its genes,
#'   background, sources, and settings and read from / written to this directory,
#'   so re-runs (e.g. across a batch) skip the network. `NULL` (the default)
#'   disables caching.
#' @param basis_gene_method How genes are assigned to factors for the primary
#'   basis-gene output (the exported `Basis_Genes_*` tables, the enrichment
#'   input, and the heatmap grouping). `"argmax"` (the default) assigns every
#'   gene to its dominant factor; `"specificity"` instead uses the sharper
#'   Kim-Park specificity markers (`NMF::extractFeatures`), leaving many genes
#'   unassigned. Changing this changes the primary scientific output.
#' @param stability_nboot Integer. If greater than 0, factor stability is
#'   assessed for each rank by refitting on `stability_nboot` random subsamples
#'   of the samples (see [nmf_stability()]); a `Stability` column is added to
#'   `consolidated_summary_df` and the per-rank tibbles are returned in
#'   `$stability`. Defaults to 0 (off), as it multiplies the fitting cost.
#' @param stability_frac Fraction of samples drawn in each stability subsample
#'   (default 0.8); only used when `stability_nboot > 0`.
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
#'     \item{`basis_genes`}{A named list of per-rank basis-gene tables (every
#'       gene to its dominant factor, by argmax).}
#'     \item{`marker_genes`}{A named list of per-rank factor-specific marker-gene
#'       tables (specificity-scored; only when `emit_marker_genes = TRUE`).}
#'     \item{`sample_assignments`}{A named list of per-rank sample-assignment
#'       tables (sample to dominant factor, with silhouette widths).}
#'     \item{`enrichment`}{A list with `per_factor`, `combined`, and `markers`
#'       named lists of per-rank g:Profiler enrichment data frames.}
#'     \item{`nmf_rds`}{A named list of file paths to the saved NMF fit objects
#'       (empty when `write_files = FALSE`).}
#'     \item{`output_dirs`}{A named list of the output directories created for
#'       the run, or `NULL` when `write_files = FALSE`.}
#'     \item{`failures`}{A data frame (columns `Rank`, `Reason`) of ranks whose
#'       NMF fit failed and were skipped; empty when every rank succeeded.}
#'     \item{`stability`}{A named list (by rank) of per-factor stability tibbles,
#'       present only when `stability_nboot > 0` (otherwise an empty list).}
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
    on_duplicate_genes = c("collapse_max", "collapse_mean", "error"),
    emit_marker_genes = FALSE,
    variance_scale = c("raw", "log"),
    custom_theme = NULL,
    factor_palette = NULL,
    umap_n_neighbors = 15,
    run_id = NULL,
    nmf_parallel = FALSE,
    nmf_cores = NULL,
    max_query_size = 10000,
    enrichment_cache = NULL,
    basis_gene_method = c("argmax", "specificity"),
    stability_nboot = 0,
    stability_frac = 0.8
) {
  on_duplicate_genes <- match.arg(on_duplicate_genes)
  variance_scale <- match.arg(variance_scale)
  basis_gene_method <- match.arg(basis_gene_method)
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
    on_duplicate_genes = on_duplicate_genes,
    variance_scale = variance_scale
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

  # --- NMF Rank Estimation Survey ---
  # Sequential execution (parallel = 0) on all platforms: NMF's parallel workers
  # cannot see the NMF namespace when it is called from another package, and it
  # is also required for reproducibility. See nmf_fit().
  cli::cli_h2("Step 2: Performing NMF Rank Estimation Survey")

  estim_real <- tryCatch({
    do.call(NMF::nmfEstimateRank, list(
      x = expr_matrix,
      range = nmf_rank,
      method = nmf_method,
      nrun = nmf_nrun,
      seed = nmf_seed,
      .options = list(parallel = 0),
      .pbackend = NA
    ))
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
  marker_genes_by_rank <- list()
  sample_assignments_by_rank <- list()
  enrichment_per_factor_by_rank <- list()
  enrichment_combined_by_rank <- list()
  marker_enrichment_by_rank <- list()
  stability_by_rank <- list()     # per-rank factor stability (P11)
  failed_ranks <- integer(0)      # ranks whose NMF fit failed (P5)
  failure_reasons <- character(0)

  for (k in nmf_rank) {

    cli::cli_rule(left = paste("Processing Rank k =", k))
    vmsg("Running NMF for k = {k} ...")
    # Capture the reason a fit fails (nmf_fit() traps the error into a warning
    # and returns NULL). The handler does not muffle the warning, so it still
    # surfaces as before; we only record its message for the failure report.
    last_fit_warning <- NA_character_
    nmf_result <- withCallingHandlers(
      run_nmf_for_rank(k, expr_matrix, nmf_method, nmf_nrun, nmf_seed, dirs$nmf_core, file_prefix,
                       write_files = write_files, nmf_parallel = nmf_parallel, n_cores = nmf_cores),
      warning = function(w) last_fit_warning <<- conditionMessage(w)
    )

    if (is.null(nmf_result)) {
      reason <- if (is.na(last_fit_warning)) "NMF fit returned NULL" else last_fit_warning
      failed_ranks <- c(failed_ranks, k)
      failure_reasons <- c(failure_reasons, reason)
      cli::cli_alert_warning("Rank k={k} failed and was skipped: {reason}")
      log_line("Rank k=", k, ": FAILED (", reason, ").")
    }

    if (!is.null(nmf_result)) {

        fits_by_rank[[as.character(k)]] <- nmf_result

        # save path (in case run_nmf_for_rank didn't)
        nmf_rds <- file.path(dirs$nmf_core, paste0("NMF_Result_Object_Rank_k", k, ".rds"))
        if (file.exists(nmf_rds)) nmf_rds_paths[[as.character(k)]] <- nmf_rds

        # Extract W and H
        W <- NMF::basis(nmf_result)
        H <- NMF::coef(nmf_result)

        # --- Assign genes to factors for the primary basis-gene output
        #     (argmax by default, or Kim-Park specificity markers). ---
        basis_genes_df <- if (identical(basis_gene_method, "specificity")) {
          nmf_marker_genes(nmf_result)
        } else {
          nmf_basis_genes(nmf_result)
        }
        basis_genes_by_rank[[as.character(k)]] <- basis_genes_df

        if (isTRUE(write_files)) {
          dir.create(dirs$basis_genes, showWarnings = FALSE, recursive = TRUE)
          readr::write_tsv(
            basis_genes_df,
            file.path(dirs$basis_genes, paste0("Basis_Genes_Rank_k", k, ".tsv"))
          )
        }

        # --- Factor stability by subsampling (opt-in, P11) ---
        if (stability_nboot > 0) {
          stab <- tryCatch(
            nmf_stability(expr_matrix, rank = k, method = nmf_method, nrun = nmf_nrun,
                          seed = nmf_seed, nboot = stability_nboot,
                          subsample_frac = stability_frac,
                          nmf_parallel = nmf_parallel, n_cores = nmf_cores),
            error = function(e) {
              cli::cli_alert_warning("Stability estimation failed for k={k}: {conditionMessage(e)}")
              NULL
            }
          )
          if (!is.null(stab)) {
            stab$Rank <- k
            stability_by_rank[[as.character(k)]] <- stab
          }
        }

        # Split into list
        basis_genes_list <- split(basis_genes_df$Gene, basis_genes_df$Factor)
        log_line("Rank k=", k, ": fitted; ", nrow(W), " genes assigned to ", k, " factors.")

        # --- Specificity-scored marker genes (additionally emitted alongside
        #     the argmax assignment; the argmax output is unchanged) ---
        if (isTRUE(emit_marker_genes)) {
          marker_genes_df <- nmf_marker_genes(nmf_result)
          marker_genes_by_rank[[as.character(k)]] <- marker_genes_df
          if (isTRUE(write_files) && nrow(marker_genes_df) > 0) {
            dir.create(dirs$basis_genes, showWarnings = FALSE, recursive = TRUE)
            readr::write_tsv(marker_genes_df, file.path(dirs$basis_genes, paste0("Marker_Genes_Rank_k", k, ".tsv")))
          }
          marker_genes_list <- split(marker_genes_df$Gene, marker_genes_df$Factor)
          marker_gost <- perform_enrichment(
            basis_genes_list = marker_genes_list,
            organism = gprofiler_organism, cutoff = gprofiler_cutoff,
            correction = gprofiler_correction, background_genes = final_nmf_genes,
            sources = gprofiler_sources, output_dir = dirs$enrichment, k = k,
            write_files = write_files, label = "Marker_Enrichment",
            max_query_size = max_query_size, enrichment_cache = enrichment_cache
          )
          marker_valid <- marker_gost[!vapply(marker_gost, is.null, logical(1))]
          marker_enrichment_by_rank[[as.character(k)]] <- dplyr::bind_rows(
            lapply(marker_valid, function(g) g$result)
          )
        }

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
          write_files = write_files,
          max_query_size = max_query_size,
          enrichment_cache = enrichment_cache
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
          write_files = write_files,
          max_query_size = max_query_size,
          enrichment_cache = enrichment_cache
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
          # Plotting is non-essential to the returned results; never let a plot
          # failure abort the whole run.
          tryCatch(
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
              nmf_seed = nmf_seed,
              user_theme = custom_theme,
              factor_palette = factor_palette,
              umap_n_neighbors = umap_n_neighbors
            ),
            error = function(e) {
              cli::cli_alert_warning("generate_rank_plots() failed for k={k}: {conditionMessage(e)}")
            }
          )
        }
      }
    }
    # --- Finalize and Global Plots ---
    cli::cli_h2("Step 4: Finalizing and Generating Global Summary Plots")
    consolidated_summary_df <- if (length(all_summaries_list) > 0) dplyr::bind_rows(all_summaries_list) else tibble::tibble()
    # Join per-factor stability onto the summary (P11), when computed.
    stability_all <- if (length(stability_by_rank) > 0) dplyr::bind_rows(stability_by_rank) else NULL
    if (!is.null(stability_all) && nrow(consolidated_summary_df) > 0) {
      consolidated_summary_df <- dplyr::left_join(consolidated_summary_df, stability_all,
                                                  by = c("Rank", "Factor"))
    }
    # Stamp the run identifier (P2) as a leading column so summaries from many
    # runs can be pooled and traced. Character, so downstream numeric-column
    # selection (e.g. the summary heatmap) ignores it.
    if (!is.null(run_id) && nrow(consolidated_summary_df) > 0) {
      consolidated_summary_df <- dplyr::mutate(consolidated_summary_df,
                                               Run_ID = as.character(run_id), .before = 1)
    }
    if (isTRUE(write_files)) {
      dir.create(dirs$summaries, showWarnings = FALSE, recursive = TRUE)
      if (nrow(consolidated_summary_df) > 0) readr::write_tsv(consolidated_summary_df, file.path(dirs$summaries, "Consolidated_Summary.tsv"))

      # generate global plots (guard if summary empty)
      tryCatch({
        generate_global_plots(all_rank_metrics_df, consolidated_summary_df, nmf_rank, nmf_nrun, dirs, file_prefix,
                              user_theme = custom_theme)
      }, error = function(e) {
        cli::cli_alert_warning("generate_global_plots() failed: {conditionMessage(e)}")
      })
    }

    time_end <- Sys.time()
    runtime <- difftime(time_end, time_start)
    cli::cli_alert_success("NMFprofileR Pipeline finished successfully in {format(runtime)}.")

    # --- Failure report (P5): ranks whose NMF fit failed and were skipped. ---
    failures <- data.frame(Rank = failed_ranks, Reason = failure_reasons,
                           stringsAsFactors = FALSE)
    if (nrow(failures) > 0) {
      cli::cli_alert_warning(
        "{nrow(failures)} of {length(nmf_rank)} rank(s) failed: {paste(failures$Rank, collapse = ', ')}."
      )
    }

    # --- Provenance: record the exact parameters, g:Profiler snapshot, and
    #     session information that produced this run (reproducibility hardening).
    run_parameters <- list(
      run_id               = run_id,
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
      marker_genes = marker_genes_by_rank,
      sample_assignments = sample_assignments_by_rank,
      enrichment = list(
        per_factor = enrichment_per_factor_by_rank,
        combined   = enrichment_combined_by_rank,
        markers    = marker_enrichment_by_rank
      ),
      nmf_rds = nmf_rds_paths,
      failures = failures,
      stability = stability_by_rank,
      output_dirs = if (isTRUE(write_files)) dirs else NULL,
      runtime = runtime,
      provenance = provenance
    )
    class(result) <- "nmf_profile"

    # --- Self-describing outputs (P2): save the whole result as a single .rds
    #     bundle and write a manifest of every file produced this run. ---
    if (isTRUE(write_files)) {
      bundle_path <- file.path(dirs$main, paste0(file_prefix, "_nmf_profile.rds"))
      manifest_tsv_path <- file.path(dirs$summaries, "manifest.tsv")
      result$provenance$files$profile_bundle <- bundle_path
      result$provenance$files$output_manifest <- manifest_tsv_path
      tryCatch({
        saveRDS(result, bundle_path)
        write_output_manifest(dirs, extra_paths = bundle_path)
      }, error = function(e) {
        cli::cli_alert_warning("Failed to write result bundle/manifest: {conditionMessage(e)}")
      })
    }

  invisible(result)
}
