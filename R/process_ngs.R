#' @title Process next generation sequencing data
#'
#' @description
#' This function performs quality control, filtering, normalization, and transformation
#' of sequencing data raw counts.
#' It can also build `phyloseq` objects for downstream ecological analyses,
#' and optionally returns intermediate processing steps.
#'
#' @details
#' - Zeros are imputed with `zCompositions::cmultRepl()` before CLR transformation.
#' - QC or other samples are removed if `remove_ids` is specified.
#' - Sample IDs in `X` and `sample_data` row names are matched and aligned.
#' - Can generate both a `phyloseq_raw` phyloseq object containing raw counts and a
#' `phyloseq_eco` object with ecosystem counts, if a `load_colname` column from
#' `sample_data` is provided to normalize the counts by microbial load (recommended best practice).
#'
#' @param X A numeric matrix or data frame of raw counts with samples as rows
#'   and features (e.g., taxa) as columns. Row names must be sample IDs.
#' @param sample_data A data frame containing sample-level data.
#'   Must include a column named `sample_id` with matching row names with `X`.
#' @param taxa_table Optional. Taxonomy annotation table to build `phyloseq` objects.
#'   Row names must match column names of `X`.
#' @param phylo_tree Optional. Phylogenetic tree to add to `phyloseq` objects.
#' @param remove_ids A regex or character vector to filter rows in `X`. Set to `NULL` to skip.
#' @param min_reads Numeric. Minimum number of total reads required per sample.
#'  Default is 500.
#' @param min_prev Numeric between 0 and 1. Minimum feature prevalence threshold.
#'  Default is 0.1 (i.e., feature must be present in >= 10 % of samples).
#' @param normalise Normalization method. One of `"load"` (microbial load data),
#'   `"TSS"` (total sum scaling), or `"none"`.
#' @param load_colname Column name in `sample_data` containing microbial load values.
#'   Required if `normalise = "load"`.
#' @param min_load Numeric. Default is 1e4. Warns if any microbial load value < min_load.
#' @param transform Transformation method. One of `"clr"` (centered log-ratio with zero imputation),
#'  `"log"` (pseudo-log using `log1p()`), or `"none"`.
#'   Note: When using `"clr"`, zero values are imputed using `zCompositions::cmultRepl()`.
#' @param impute_control A named list of arguments to be passed to `zCompositions::cmultRepl()`.
#' @param return_all Logical. If `TRUE`, additional intermediate data matrices
#'   (`X_matched`, `X_norm`, `X_prev`) are included in the output. Default is `FALSE`.
#' @param raw_phyloseq Logical. If `TRUE`, constructs a `phyloseq` object with
#'   the table of raw counts (filtered failed runs if needed). Default is `TRUE`.
#' @param eco_phyloseq Logical. If `TRUE`, constructs a `phyloseq` object
#'   with the ecosystem abundances (i.e. after `normalise = "load"`). Default is `TRUE`.
#' @param verbose Logical. If `TRUE`, prints progress messages during execution.
#'   Default is `TRUE`.
#'
#' @return A named list containing:
#' \describe{
#'   \item{`X_processed`}{Matrix of processed feature counts after filtering,
#'                        normalization, and transformation.}
#'   \item{`sdata_final`}{Matched and filtered `sample_data` corresponding to
#'                        retained samples.}
#'   \item{`phyloseq_raw`}{`phyloseq` object created from raw filtered data.
#'                         `NULL` if `raw_phyloseq = FALSE`.}
#'   \item{`phyloseq_eco`}{`phyloseq` object from ecosystem abundance data.
#'                         `NULL` if `eco_phyloseq = FALSE` or `normalise != "load"`.}
#'   \item{`X_matched`}{(Optional) Matched and filtered count matrix, pre-normalization.
#'                      Returned only if `return_all = TRUE`.}
#'   \item{`X_norm`}{(Optional) Normalized count matrix. Returned only if `return_all = TRUE`.}
#'   \item{`X_prev`}{(Optional) Prevalence-filtered matrix, pre-transformation.
#'                   Returned only if `return_all = TRUE`.}
#' }
#'
#' @seealso
#' * [build_phyloseq()]
#' * [zCompositions::cmultRepl()]
#'
#' @references
#' #' McMurdie, P. J., & Holmes, S. (2013).
#' phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data.
#' *PLoS ONE*, 8(4), e61217. \doi{10.1371/journal.pone.0061217}
#'
#' Martín-Fernández, J. A., Hron, K., Templ, M., Filzmoser, P., & Palarea-Albaladejo, J. (2015).
#' Bayesian-multiplicative treatment of count zeros in compositional data sets.
#' *Statistical Modelling*, 15(2), 134–158. \doi{10.1177/1471082X14535524}
#'
#' Palarea-Albaladejo, J., & Martín-Fernández, J. A. (2015).
#' zCompositions—R package for multivariate imputation of left-censored data under a compositional approach.
#' *Chemometrics and Intelligent Laboratory Systems*, 143, 85–96. \doi{10.1016/j.chemolab.2015.02.019}
#'
#' Gloor, G. B., Macklaim, J. M., Pawlowsky-Glahn, V., & Egozcue, J. J. (2017).
#' Microbiome datasets are compositional: And this is not optional.
#' *Frontiers in Microbiology*, 8, 2224. \doi{10.3389/fmicb.2017.02224}
#'
#' Vandeputte, D., Kathagen, G., D’hoe, K., Vieira-Silva, S., Valles-Colomer, M., Sabino, J.,
#' Wang, J., Tito, R. Y., De Commer, L., Darzi, Y., Vermeire, S., Falony, G., & Raes, J. (2017).
#' Quantitative microbiome profiling links gut community variation to microbial load. *Nature*,
#' 551(7681), 507–511. \doi{10.1038/nature24460}
#'
#' @examples
#' if (requireNamespace("phyloseq", quietly = TRUE)) {
#' mock_X <- matrix(sample(0:1000, 25, replace = TRUE),
#'                  nrow = 5,
#'                  dimnames = list(paste0("sample", 1:5),
#'                  paste0("ASV", 1:5))
#'                  )
#'
#' mock_sample_data <- data.frame(
#'   sample_id = paste0("sample", 1:5),
#'   load = c(1e5, 2e5, 1e4, 5e4, 1.5e5),
#'   condition = factor(rep(c("A", "B"), length.out = 5)),
#'   row.names = paste0("sample", 1:5)
#'   )
#'
#' mock_taxa_table <- data.frame(
#'   Kingdom = rep("Bacteria", 5),
#'   Genus = paste0("Genus", 1:5),
#'   row.names = paste0("ASV", 1:5)
#'   )
#'
#' result <- process_ngs(
#'   X = mock_X,
#'   sample_data = mock_sample_data,
#'   taxa_table = mock_taxa_table,
#'   normalise = "load",
#'   load_colname = "load",
#'   transform = "none",
#'   verbose = FALSE
#'   )
#' }
#'
#' @export
process_ngs <- function(X, sample_data, taxa_table = NULL,
                        phylo_tree = NULL, remove_ids = NULL, min_reads = 500,
                        min_prev = 0.1, normalise = c("load", "TSS", "none"),
                        load_colname = NULL, min_load = 1e4,
                        transform = c("clr", "log", "none"),
                        impute_control = list(method = "GBM",
                                              output = "p-counts",
                                              z.delete = FALSE,
                                              z.warning = 1,
                                              suppress.print = TRUE),
                        raw_phyloseq = TRUE, eco_phyloseq = TRUE,
                        return_all = FALSE, verbose = TRUE) {
  # Match or set arguments
  normalise <- match.arg(normalise)
  transform <- match.arg(transform)

  phyloseq_raw <- NULL #updated if raw_phyloseq = TRUE
  phyloseq_eco <- NULL #unless data normalised by microbial load

  # Check X
  X <- check_omics_matrix(X, platform = "ngs")

  # Check column and row names
  check_ids <- check_sample_ids(X, sample_data, do_match = FALSE,
                                verbose = FALSE)
  X <- check_ids$X
  sample_data <- check_ids$sample_data

  # Check parameters
  if (normalise == "load") {
    if (is.null(load_colname) || !is.character(load_colname) || length(load_colname) != 1) {
      stop("'load_colname' must be a single character string specifying a column in 'sample_data'
           containing microbial load values.")
    }
    if (!load_colname %in% colnames(sample_data)) {
      stop(paste0("Column '", load_colname, "' not found in 'sample_data'.\n",
                  "Please provide a valid column name for microbial load"))
    }
  }

  if (eco_phyloseq && normalise != "load") {
    stop("Phyloseq object with ecosystem abundances cannot be built unless
         normalise is set to 'load'")
  }

  # Phyloseq objects
  if (any(raw_phyloseq, eco_phyloseq)) {
    check_ngs <- check_ngs_input(X, taxa_table, phylo_tree, verbose = verbose)
    X <- check_ngs$X
    taxa_table <- check_ngs$taxa_table
  }

  # 1. Remove samples (e.g. QCs, blanks, etc.)
  X_noqc <- remove_samples(X, remove_ids = remove_ids)

  # 2. Filter failed samples (very low sequencing depth)
  n_reads <- rowSums(X_noqc)
  failed_samples <- which(n_reads == 0 | n_reads < min_reads)
  if (length(failed_samples) != 0) {
    X_filtered <- X_noqc[-failed_samples, ]
    if (verbose) {
      message(length(failed_samples),
              " samples removed due to zero total read count or < ",
              min_reads, " 'min_reads'.\n")
    }
  } else {
    X_filtered <- X_noqc
  }
  if (nrow(X_filtered) == 0) {
    stop("No samples remain after row names and min_reads filtering.")
  }

  # Match rownames across tables
  match_ids <- check_sample_ids(X_filtered, sample_data, do_match = TRUE,
                                verbose = verbose)
  X_matched <- match_ids$X
  sdata_matched <- match_ids$sample_data

  # Build phyloseq of the raw matched counts (needed for alpha diversity measures)
  if (raw_phyloseq) {
    phyloseq_raw <- build_phyloseq(X_matched, sdata_matched,
                                   taxa_table, phylo_tree,
                                   taxa_in_rows= FALSE, verbose = FALSE)
  }

  # 3. Normalisation
  if (normalise == "none") {
    X_norm <- X_matched
  } else {
    X_tss <- sweep(X_matched, 1, rowSums(X_matched), FUN = "/")
    if (normalise == "TSS") {
      X_norm <- X_tss
    } else if (normalise == "load") {
      microbial_load <- sdata_matched[[load_colname]]
      if (any(is.na(microbial_load), microbial_load == 0)) {
        warning("Microbial load column contains NA or zeros.\n",
                "Inspect microbial load or remove samples where load could not be assessed.\n",
                immediate. = TRUE)
      } else if (any(microbial_load < min_load)) {
        warning("Microbial load for some samples is <", min_load," 'min_load'.\n",
                "Double-check load values unless this is expected.\n",
                immediate. = TRUE)
      }
      X_norm <- sweep(X_tss, 1, microbial_load, FUN = "*")
      if (eco_phyloseq) {
        phyloseq_eco <- build_phyloseq(X_norm, sdata_matched,
                                       taxa_table, phylo_tree,
                                       taxa_in_rows= FALSE, verbose = FALSE)
      }
    }
  }

  # 4. Prevalence filtering
  X_prev <- remove_features(X_norm, label = "0", min_prev= min_prev,
                            verbose = verbose)

  # 5. Transformation
  if (transform == "none") {
    X_trans <- X_prev
  } else if (transform == "clr") {
    if (verbose) message("Zeros will be imputed with 'zCompositions::cmultRepl()' prior to clr-transformation.\n")

    imp_args <- utils::modifyList(impute_control, list(X = X_prev))

    X_imp <- do.call(zCompositions::cmultRepl, imp_args)

    X_trans <- t(apply(X_imp, 1, function(x) {log(x) - mean(log(x))}))
  } else if (transform == "log") {
    if (verbose) message("Pseudo-log transform log1p() will be used for log-transformation.\n")
    X_trans <- log1p(X_prev)
  }

  output <- list(X_processed = X_trans,
                 sdata_final = sdata_matched,
                 phyloseq_raw = phyloseq_raw,
                 phyloseq_eco = phyloseq_eco)

  if (return_all) {
    output$X_matched <- X_matched
    output$X_norm <- X_norm
    output$X_prev <- X_prev
  }

  return(output)
}
