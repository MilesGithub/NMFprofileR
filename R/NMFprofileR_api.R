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
#' @param nmf_parallel Logical; if `FALSE` (the default) the `nrun` consensus
#'   runs are executed sequentially. If `TRUE` they are spread across a local
#'   PSOCK cluster. The parallel path is opt-in because NMF's default parallel
#'   backend fails when NMF is invoked from inside another package (its workers
#'   cannot see the NMF namespace); here the cluster is configured to load NMF on
#'   each worker, which both fixes that and, because NMF pre-generates the RNG
#'   stream for every run, yields results identical to the sequential path for a
#'   fixed `seed`. Falls back to sequential (with a warning) if the cluster
#'   cannot be created.
#' @param n_cores Number of worker processes to use when `nmf_parallel = TRUE`.
#'   `NULL` (the default) uses one fewer than the number of detected cores,
#'   capped at `nrun`.
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
nmf_fit <- function(expr_matrix, rank, method = "brunet", nrun = 20, seed = 123456,
                    nmf_parallel = FALSE, n_cores = NULL) {
  # Sequential fit: parallel = 0, no backend. This is the safe default -- NMF's
  # parallel execution otherwise requires the NMF package to be attached in the
  # worker processes, which is not the case when NMF is called from another
  # package's namespace (it fails there with "none of the packages are loaded").
  run_sequential <- function() {
    NMF::nmf(
      x = expr_matrix, rank = rank, method = method, nrun = nrun, seed = seed,
      .options = list(parallel = 0), .pbackend = NA
    )
  }

  # Parallel fit: NMF runs the consensus runs across its own local cluster when
  # given `.options = list(parallel = n)`. It selects the packages to load on
  # the workers from the ATTACHED search path, so when NMFprofileR calls this
  # from its own namespace (NMF loaded but not attached) NMF must be attached for
  # the duration of the fit -- otherwise the workers fail with "none of the
  # packages are loaded". We attach it via attachNamespace() (not library/require
  # so it is valid in package code), only if needed, and restore the search path
  # afterwards. NMF seeds each run deterministically (via rngtools), so for a
  # fixed seed the parallel result is identical to the sequential path.
  run_parallel <- function() {
    n <- if (is.null(n_cores)) max(1L, parallel::detectCores() - 1L) else as.integer(n_cores)
    n <- min(n, nrun)
    if (is.na(n) || n < 2L) return(run_sequential())

    if (!"package:NMF" %in% search()) {
      tryCatch({
        suppressPackageStartupMessages(attachNamespace("NMF"))
        on.exit(try(detach("package:NMF", character.only = TRUE), silent = TRUE), add = TRUE)
      }, error = function(e) NULL)
    }

    NMF::nmf(
      x = expr_matrix, rank = rank, method = method, nrun = nrun, seed = seed,
      .options = list(parallel = n)
    )
  }

  tryCatch({
    if (isTRUE(nmf_parallel) && nrun > 1L) {
      tryCatch(
        run_parallel(),
        error = function(e) {
          warning(sprintf("Parallel NMF failed (%s); falling back to sequential.",
                          conditionMessage(e)), call. = FALSE)
          run_sequential()
        }
      )
    } else {
      run_sequential()
    }
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
#' @param max_query_size Integer; gene sets larger than this are skipped rather
#'   than sent to g:Profiler (default 10000).
#' @param enrichment_cache Optional path to a directory used to cache g:Profiler
#'   results by query hash; `NULL` (the default) disables caching. See
#'   [NMFprofileR()].
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
                           cutoff = 0.05,
                           max_query_size = 10000,
                           enrichment_cache = NULL) {
  per_factor <- perform_enrichment(
    basis_genes_list = basis_genes_list, organism = organism, cutoff = cutoff,
    correction = correction, background_genes = background_genes, sources = sources,
    output_dir = NULL, k = NA, write_files = FALSE,
    max_query_size = max_query_size, enrichment_cache = enrichment_cache
  )
  combined <- perform_combined_enrichment(
    basis_genes_list = basis_genes_list, organism = organism, cutoff = cutoff,
    correction = correction, background_genes = background_genes, sources = sources,
    output_dir = NULL, k = NA, write_files = FALSE,
    max_query_size = max_query_size, enrichment_cache = enrichment_cache
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

#' Extract a basis matrix (genes x factors) from a fit or matrix
#'
#' @param x An `NMFfit` object or a numeric basis matrix.
#' @return A numeric matrix with gene rownames and factor colnames.
#' @noRd
as_basis_matrix <- function(x) {
  # Accept a plain matrix, or any NMF fit object (whose concrete S4 class varies,
  # e.g. NMFfitX1); NMF::basis() extracts the basis matrix from the latter.
  W <- if (is.matrix(x)) x else tryCatch(NMF::basis(x), error = function(e) as.matrix(x))
  if (is.null(colnames(W))) colnames(W) <- paste0("Factor_", seq_len(ncol(W)))
  W
}

#' Align NMF factors across runs by basis-loading correlation
#'
#' Matches the factors of several NMF runs to a common reference so that "the
#' same" factor can be tracked across runs (e.g. across cohorts, ranks, or
#' resamples). The first element of `x` is the reference; every other run's
#' factors are correlated against the reference factors on their shared genes and
#' matched one-to-one.
#'
#' @param x A named list (length >= 2) of `NMFfit` objects and/or basis matrices
#'   (genes in rows with gene-symbol rownames, factors in columns). The names
#'   identify the runs; the first element is the reference.
#' @param method How to pick the one-to-one matches: `"greedy"` (the default)
#'   repeatedly takes the highest remaining correlation; `"hungarian"` finds the
#'   globally optimal assignment via `clue::solve_LSAP()` (used only if the
#'   `clue` package is installed, otherwise it falls back to greedy).
#' @param min_cor Minimum absolute correlation for a pair to be reported as a
#'   match (default 0.5).
#' @param cor_method Correlation method, `"pearson"` (default) or `"spearman"`.
#'
#' @return A list with:
#'   \describe{
#'     \item{`reference`}{The name of the reference run.}
#'     \item{`correlations`}{A tibble (`run`, `factor`, `ref_factor`,
#'       `correlation`) of every run factor against every reference factor.}
#'     \item{`matches`}{A tibble (`run`, `factor`, `ref_factor`, `correlation`)
#'       of the selected one-to-one matches with `abs(correlation) >= min_cor`.}
#'   }
#' @seealso [nmf_fit()]
#' @export
#'
#' @examples
#' \dontrun{
#' fit_a <- nmf_fit(m, rank = 3)
#' fit_b <- nmf_fit(m, rank = 3, seed = 99)
#' aln <- align_factors(list(run_a = fit_a, run_b = fit_b))
#' aln$matches
#' }
align_factors <- function(x,
                          method = c("greedy", "hungarian"),
                          min_cor = 0.5,
                          cor_method = c("pearson", "spearman")) {
  method <- match.arg(method)
  cor_method <- match.arg(cor_method)

  if (!is.list(x) || length(x) < 2) {
    stop("`x` must be a list of at least two runs (NMFfit objects or basis matrices).",
         call. = FALSE)
  }
  run_names <- names(x)
  if (is.null(run_names) || any(!nzchar(run_names)) || anyDuplicated(run_names)) {
    stop("`x` must have unique, non-empty names identifying the runs.", call. = FALSE)
  }

  bases <- lapply(x, as_basis_matrix)
  ref_name <- run_names[1]
  ref <- bases[[1]]

  greedy_match <- function(cormat) {
    # cormat: query factors (rows) x reference factors (cols)
    pairs <- list()
    remaining_rows <- rownames(cormat)
    remaining_cols <- colnames(cormat)
    while (length(remaining_rows) > 0 && length(remaining_cols) > 0) {
      sub <- cormat[remaining_rows, remaining_cols, drop = FALSE]
      idx <- which.max(abs(sub))
      rc <- arrayInd(idx, dim(sub))
      r <- remaining_rows[rc[1]]; cc <- remaining_cols[rc[2]]
      pairs[[length(pairs) + 1]] <- data.frame(
        factor = r, ref_factor = cc, correlation = sub[rc[1], rc[2]],
        stringsAsFactors = FALSE
      )
      remaining_rows <- setdiff(remaining_rows, r)
      remaining_cols <- setdiff(remaining_cols, cc)
    }
    do.call(rbind, pairs)
  }

  hungarian_match <- function(cormat) {
    # Maximize total abs correlation -> minimize its negative. solve_LSAP needs a
    # square-ish nonnegative cost; pad to square with a high cost.
    cost <- max(abs(cormat)) - abs(cormat)
    nr <- nrow(cost); nc <- ncol(cost)
    n <- max(nr, nc)
    padded <- matrix(max(cost) + 1, n, n)
    padded[seq_len(nr), seq_len(nc)] <- cost
    assign <- clue::solve_LSAP(padded)
    out <- lapply(seq_len(nr), function(i) {
      j <- assign[i]
      if (j > nc) return(NULL)
      data.frame(factor = rownames(cormat)[i], ref_factor = colnames(cormat)[j],
                 correlation = cormat[i, j], stringsAsFactors = FALSE)
    })
    do.call(rbind, out)
  }

  use_hungarian <- identical(method, "hungarian") && requireNamespace("clue", quietly = TRUE)

  cor_long <- list()
  match_long <- list()

  for (i in 2:length(bases)) {
    run <- run_names[i]
    q <- bases[[i]]
    shared <- intersect(rownames(ref), rownames(q))
    if (length(shared) < 2) {
      cli::cli_alert_warning("Run {run}: fewer than 2 shared genes with the reference; skipping.")
      next
    }
    cormat <- stats::cor(q[shared, , drop = FALSE], ref[shared, , drop = FALSE],
                         method = cor_method)
    # long correlations
    cl <- expand.grid(factor = rownames(cormat), ref_factor = colnames(cormat),
                      stringsAsFactors = FALSE)
    cl$run <- run
    cl$correlation <- mapply(function(a, b) cormat[a, b], cl$factor, cl$ref_factor)
    cor_long[[run]] <- cl[, c("run", "factor", "ref_factor", "correlation")]

    m <- if (use_hungarian) hungarian_match(cormat) else greedy_match(cormat)
    if (!is.null(m) && nrow(m) > 0) {
      m <- m[abs(m$correlation) >= min_cor, , drop = FALSE]
      if (nrow(m) > 0) {
        m$run <- run
        match_long[[run]] <- m[, c("run", "factor", "ref_factor", "correlation")]
      }
    }
  }

  empty <- tibble::tibble(run = character(), factor = character(),
                          ref_factor = character(), correlation = numeric())
  list(
    reference = ref_name,
    correlations = if (length(cor_long)) tibble::as_tibble(do.call(rbind, cor_long)) else empty,
    matches = if (length(match_long)) tibble::as_tibble(do.call(rbind, match_long)) else empty
  )
}

#' Project new samples onto learned NMF factors
#'
#' Scores new samples against an existing factorization by solving for their
#' coefficient matrix H with the basis matrix W held fixed
#' (`min ||V - W H||` s.t. `H >= 0`), using non-negative multiplicative updates.
#' Genes are aligned on the shared symbols between W and the new data.
#'
#' @param fit An `NMFfit` object (or a numeric basis matrix W, genes x factors)
#'   defining the learned factors.
#' @param newdata A numeric matrix or data frame of new samples (genes in rows
#'   with gene-symbol rownames, samples in columns).
#' @param n_iter Maximum number of multiplicative-update iterations (default 200).
#' @param tol Convergence tolerance on the maximum change in H (default 1e-6).
#'
#' @return A data frame with one row per new sample: `SampleID`, one coefficient
#'   column per factor (`Factor_1 ... Factor_k`), and `Dominant_Factor`.
#' @seealso [nmf_fit()], [nmf_sample_assignments()]
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- nmf_fit(m, rank = 3)
#' nmf_project(fit, newdata = m_new)
#' }
nmf_project <- function(fit, newdata, n_iter = 200L, tol = 1e-6) {
  W <- as_basis_matrix(fit)
  V <- as.matrix(newdata)
  if (is.null(rownames(V))) stop("`newdata` must have gene symbols as rownames.", call. = FALSE)

  shared <- intersect(rownames(W), rownames(V))
  if (length(shared) < 2) {
    stop("Fewer than 2 shared genes between the fit and `newdata`.", call. = FALSE)
  }
  n_missing <- nrow(W) - length(shared)
  if (n_missing > 0) {
    cli::cli_alert_warning("{n_missing} of {nrow(W)} basis genes are absent from `newdata`; projecting on {length(shared)} shared genes.")
  }

  Ws <- W[shared, , drop = FALSE]
  Vs <- V[shared, , drop = FALSE]
  Vs[!is.finite(Vs)] <- 0
  Vs[Vs < 0] <- 0

  k <- ncol(Ws)
  n <- ncol(Vs)
  # Deterministic non-negative initialisation (no RNG) so projection is reproducible.
  H <- matrix(1, nrow = k, ncol = n)
  WtW <- crossprod(Ws)          # k x k
  WtV <- crossprod(Ws, Vs)      # k x n
  eps <- .Machine$double.eps
  for (i in seq_len(n_iter)) {
    H_new <- H * (WtV / (WtW %*% H + eps))
    if (max(abs(H_new - H)) < tol) { H <- H_new; break }
    H <- H_new
  }

  factor_names <- if (is.null(colnames(Ws))) paste0("Factor_", seq_len(k)) else colnames(Ws)
  score_mat <- t(H)
  colnames(score_mat) <- factor_names
  sample_ids <- if (is.null(colnames(Vs))) paste0("Sample_", seq_len(n)) else colnames(Vs)

  dominant <- factor_names[max.col(score_mat, ties.method = "first")]
  scores <- as.data.frame(score_mat)
  out <- data.frame(SampleID = sample_ids, scores, Dominant_Factor = dominant,
                    stringsAsFactors = FALSE, check.names = FALSE)
  rownames(out) <- NULL
  out
}

#' Assess NMF factor stability by subsampling
#'
#' Quantifies how reproducible each factor is by refitting NMF on random subsets
#' of the samples and measuring how well the resampled factors recover the
#' reference factors (fitted on the full matrix). For each resample the basis
#' loadings are correlated against the reference and each reference factor's best
#' match is recorded; the per-factor stability is the mean of these best matches
#' across resamples (1 = perfectly reproducible).
#'
#' @param expr_matrix A preprocessed, non-negative numeric matrix (genes x
#'   samples), as from [nmf_preprocess()].
#' @param rank The factorization rank (integer).
#' @param method,nrun,seed Passed to [nmf_fit()] for every fit.
#' @param nboot Number of subsampling iterations (default 30).
#' @param subsample_frac Fraction of samples drawn (without replacement) in each
#'   iteration (default 0.8).
#' @param nmf_parallel,n_cores Passed to [nmf_fit()].
#'
#' @return A tibble with columns `Factor` and `Stability` (mean best-match
#'   correlation across resamples), one row per reference factor.
#' @seealso [nmf_fit()], [nmf_rank_diagnostics()]
#' @export
#'
#' @examples
#' \dontrun{
#' m <- nmf_preprocess(example_expression_data, expression_threshold = 10,
#'                     variance_quantile = 0.9)
#' nmf_stability(m, rank = 3, nrun = 10, nboot = 20)
#' }
nmf_stability <- function(expr_matrix, rank, method = "brunet", nrun = 10, seed = 123456,
                          nboot = 30, subsample_frac = 0.8,
                          nmf_parallel = FALSE, n_cores = NULL) {
  factor_names <- paste0("Factor_", seq_len(rank))
  na_result <- tibble::tibble(Factor = factor_names, Stability = NA_real_)

  ref_fit <- nmf_fit(expr_matrix, rank = rank, method = method, nrun = nrun, seed = seed,
                     nmf_parallel = nmf_parallel, n_cores = n_cores)
  if (is.null(ref_fit)) {
    cli::cli_alert_warning("Reference fit failed; returning NA stability.")
    return(na_result)
  }
  W_ref <- NMF::basis(ref_fit)

  n_samples <- ncol(expr_matrix)
  size <- max(rank + 1L, floor(subsample_frac * n_samples))
  if (size >= n_samples) {
    cli::cli_alert_warning("subsample_frac leaves no samples out; stability is uninformative.")
  }

  best_per_boot <- vector("list", nboot)
  for (b in seq_len(nboot)) {
    set.seed(seed + b)
    cols <- sample.int(n_samples, size = min(size, n_samples))
    fit_b <- nmf_fit(expr_matrix[, cols, drop = FALSE], rank = rank, method = method,
                     nrun = nrun, seed = seed, nmf_parallel = nmf_parallel, n_cores = n_cores)
    if (is.null(fit_b)) next
    W_b <- NMF::basis(fit_b)
    shared <- intersect(rownames(W_ref), rownames(W_b))
    if (length(shared) < 2) next
    cormat <- abs(stats::cor(W_ref[shared, , drop = FALSE], W_b[shared, , drop = FALSE]))
    best_per_boot[[b]] <- apply(cormat, 1, max)   # best match per reference factor
  }

  best_per_boot <- best_per_boot[!vapply(best_per_boot, is.null, logical(1))]
  if (length(best_per_boot) == 0) return(na_result)

  stability <- rowMeans(do.call(cbind, best_per_boot))
  tibble::tibble(Factor = factor_names, Stability = as.numeric(stability))
}

#' Run NMFprofileR over many cohorts and consolidate the results
#'
#' A batch driver that applies [NMFprofileR()] to each of several expression
#' matrices ("cohorts"), writing each cohort's outputs to its own subdirectory of
#' `output_dir` and returning a single consolidated per-factor summary across all
#' cohorts. It is resilient (one cohort failing does not abort the batch),
#' resumable (already-completed cohorts can be skipped), and tags every row with
#' the cohort's `run_id` so results can be pooled and traced.
#'
#' @param cohorts A named list of expression matrices or data frames (one per
#'   cohort), in the form accepted by [NMFprofileR()]. The names are used both as
#'   the `run_id` and as the output subdirectory for each cohort.
#' @param output_dir A directory under which each cohort's `<name>_Results` tree
#'   is written.
#' @param ... Further arguments passed to [NMFprofileR()] (e.g. `nmf_rank`,
#'   `nmf_nrun`, `gprofiler_sources`, `nmf_parallel`, `enrichment_cache`). Do not
#'   pass `expression_data`, `output_prefix`, or `run_id`; the driver sets those
#'   per cohort.
#' @param skip_existing Logical. If `TRUE` (the default), a cohort whose results
#'   directory already exists is not re-run; its summary is loaded from disk and
#'   still included in the consolidated output.
#' @param on_error One of `"continue"` (the default; log the failing cohort and
#'   carry on) or `"stop"` (abort the whole batch on the first cohort error).
#' @param write_batch_summary Logical. If `TRUE` (the default) the consolidated
#'   summary is also written to `output_dir/Batch_Consolidated_Summary.tsv`.
#'
#' @return Invisibly, a list with:
#'   \describe{
#'     \item{`consolidated`}{A data frame binding every cohort's
#'       `consolidated_summary_df` (each tagged with its `Run_ID`).}
#'     \item{`failures`}{A data frame (`Cohort`, `Reason`) of cohorts that errored.}
#'     \item{`skipped`}{A character vector of cohorts skipped because their
#'       results already existed.}
#'   }
#' @seealso [NMFprofileR()]
#' @export
#'
#' @examples
#' \dontrun{
#' data("example_expression_data")
#' # Split the samples into two illustrative cohorts.
#' half <- ncol(example_expression_data) %/% 2
#' cohorts <- list(
#'   cohort_A = example_expression_data[, seq_len(half)],
#'   cohort_B = example_expression_data[, (half + 1):ncol(example_expression_data)]
#' )
#' batch <- run_nmf_batch(cohorts, output_dir = "batch_out", nmf_rank = 2:3, nmf_nrun = 10)
#' batch$consolidated
#' }
run_nmf_batch <- function(cohorts,
                          output_dir,
                          ...,
                          skip_existing = TRUE,
                          on_error = c("continue", "stop"),
                          write_batch_summary = TRUE) {
  on_error <- match.arg(on_error)

  if (!is.list(cohorts) || length(cohorts) == 0) {
    stop("`cohorts` must be a non-empty named list of expression matrices/data frames.",
         call. = FALSE)
  }
  nms <- names(cohorts)
  if (is.null(nms) || any(is.na(nms)) || any(!nzchar(nms)) || anyDuplicated(nms)) {
    stop("`cohorts` must have unique, non-empty names (used as run ids and output dirs).",
         call. = FALSE)
  }
  dots <- list(...)
  reserved <- intersect(names(dots), c("expression_data", "output_prefix", "run_id"))
  if (length(reserved) > 0) {
    stop("Do not pass ", paste(sprintf("`%s`", reserved), collapse = ", "),
         " to run_nmf_batch(); the driver sets these per cohort.", call. = FALSE)
  }

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Load a completed cohort's consolidated summary from disk (bundle preferred).
  load_cohort_summary <- function(results_dir, run_id) {
    bundle <- list.files(results_dir, pattern = "_nmf_profile\\.rds$", full.names = TRUE)
    if (length(bundle) >= 1) {
      prof <- tryCatch(readRDS(bundle[[1]]), error = function(e) NULL)
      if (!is.null(prof) && !is.null(prof$consolidated_summary_df)) {
        return(prof$consolidated_summary_df)
      }
    }
    csv <- file.path(results_dir, "Summaries", "Consolidated_Summary.tsv")
    if (file.exists(csv)) {
      df <- tryCatch(readr::read_tsv(csv, show_col_types = FALSE), error = function(e) NULL)
      if (!is.null(df)) {
        if (!"Run_ID" %in% names(df)) df <- dplyr::mutate(df, Run_ID = run_id, .before = 1)
        return(df)
      }
    }
    NULL
  }

  cli::cli_h1("Running NMF batch over {length(cohorts)} cohort(s)")

  summaries <- list()
  skipped <- character(0)
  failed_cohorts <- character(0)
  failure_reasons <- character(0)

  for (i in seq_along(cohorts)) {
    name <- nms[i]
    safe_name <- gsub("[^A-Za-z0-9._-]+", "_", name)
    output_prefix <- file.path(output_dir, safe_name)
    results_dir <- paste0(output_prefix, "_Results")

    cli::cli_rule(left = "Cohort {i}/{length(cohorts)}: {name}")

    if (isTRUE(skip_existing) && dir.exists(results_dir)) {
      cli::cli_alert_info("Skipping {name}: results already exist at {.path {results_dir}}.")
      loaded <- load_cohort_summary(results_dir, name)
      if (!is.null(loaded)) summaries[[name]] <- loaded
      skipped <- c(skipped, name)
      next
    }

    res <- tryCatch(
      do.call(NMFprofileR, c(
        list(expression_data = cohorts[[name]], output_prefix = output_prefix, run_id = name),
        dots
      )),
      error = function(e) e
    )

    if (inherits(res, "error")) {
      msg <- conditionMessage(res)
      if (identical(on_error, "stop")) {
        stop(sprintf("Cohort '%s' failed: %s", name, msg), call. = FALSE)
      }
      cli::cli_alert_danger("Cohort {name} failed: {msg}")
      failed_cohorts <- c(failed_cohorts, name)
      failure_reasons <- c(failure_reasons, msg)
      next
    }

    summaries[[name]] <- res$consolidated_summary_df
  }

  consolidated <- if (length(summaries) > 0) dplyr::bind_rows(summaries) else tibble::tibble()
  failures <- data.frame(Cohort = failed_cohorts, Reason = failure_reasons,
                         stringsAsFactors = FALSE)

  if (isTRUE(write_batch_summary) && nrow(consolidated) > 0) {
    out_path <- file.path(output_dir, "Batch_Consolidated_Summary.tsv")
    tryCatch(
      readr::write_tsv(consolidated, out_path),
      error = function(e) cli::cli_alert_warning("Failed to write batch summary: {conditionMessage(e)}")
    )
  }

  cli::cli_alert_success(
    "Batch complete: {length(summaries)} cohort(s) summarized ({length(skipped)} skipped, {nrow(failures)} failed)."
  )

  invisible(list(consolidated = consolidated, failures = failures, skipped = skipped))
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
  if (!is.null(x$failures) && nrow(x$failures) > 0) {
    cat(sprintf("  Ranks failed : %s\n", paste(x$failures$Rank, collapse = ", ")))
  }
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