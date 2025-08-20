#' @title Validate and prepare an omics data matrix
#'
#' @description
#' This internal helper function checks that the input matrix `X` is numeric,
#' does not contain NaN values, and, depending on the platform, does not contain negative values.
#' It also verifies that row names (typically sample IDs) are present.
#'
#' @param X A numeric matrix or data frame representing omics data.
#'  Samples should be in rows and features in columns.
#' @param platform Character string specifying the omics platform.
#'  Must be one of `"ms"` (mass spectrometry), `"nmr"` (nuclear magnetic resonance),
#'  or `"ngs"` (next-generation sequencing). Negative values are only allowed if `platform = "nmr"`.
#'
#' @return A numeric matrix with validated structure and values.
#'
#' @keywords internal
#'
#' @noRd
check_omics_matrix <- function(X, platform = c("ms", "nmr", "ngs")) {
  platform <- match.arg(platform)
  X <- as.matrix(X)
  if (!is.numeric(X)) {
    stop("'X' must be a numeric matrix or data frame.")
  }
  if (any(is.nan(X))) {
    stop("'X' contains NaN values.")
  }
  if (any(X < 0, na.rm = TRUE) && platform != "nmr") {
    stop("'X' contains negative values.")
  }
  if (is.null(rownames(X))) {
    stop("No `X` row names found: 'sample_id' from 'sample_data' needs to be set as row names.")
  }
  return(X)
}

#' @title Match and align sample identifiers between omics matrix and sample-level data tables
#'
#' @description
#' This internal helper function checks and aligns sample identifiers between
#' an omics data matrix (`X`) and the corresponding sample data (`sample_data`).
#' It ensures consistent and properly set `sample_id`, and optionally aligns and
#' subsets both tables to shared samples.
#'
#' @param X A numeric matrix or data frame with sample IDs as row names.
#' @param sample_data A data frame containing sample data. Must include a column
#'  named `"sample_id"` which matches its row names.
#' @param do_match Logical. If `TRUE`, both `X` and `sample_data` are subset and
#'  aligned to shared sample IDs, following `sample_data` row order.
#' @param verbose Logical. If `TRUE`, prints a message indicating the number of matched samples.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{`X`}{Omics matrix with matched rownames if `do_match = TRUE`.}
#'   \item{`sample_data`}{Sample data with matched rownames if `do_match = TRUE`.}
#'
#' }
#'
#' @keywords internal
#'
#' @noRd
check_sample_ids <- function(X, sample_data, do_match = TRUE,
                             verbose = TRUE) {
  # Check labels in sample_data
  if (!"sample_id" %in% colnames(sample_data)) {
    stop("'sample_data' must include a 'sample_id' column of unique labels used as row names.")
  } else {
    sample_ids <- as.character(sample_data$sample_id)
    if (length(unique(sample_ids)) != nrow(sample_data)) {
      stop("Ensure 'sample_id' has", nrow(sample_data),
           " unique labels and is set as row names in 'X' and 'sample_data'.")
    } else if (!identical(sample_ids, rownames(sample_data))) {
      warning("'sample_id' column was found in 'sample_data' and will be set as row names.\n",
              immediate. = TRUE)
      rownames(sample_data) <- sample_ids
    }
  }

  # Check ids overlap
  common_ids <- intersect(rownames(X), rownames(sample_data))
  if (length(common_ids) == 0) {
    stop("Ensure 'sample_id' column from 'sample_data' is set as row names in 'X'.")
  }

  if (verbose) {
    message(length(common_ids), " samples found in common between ",
            nrow(X), " rows in 'X' and ",
            nrow(sample_data), " rows in 'sample_data'.\n")
  }

  # Generate matched tables
  if (do_match) {
    sdata_matched <- droplevels(sample_data[rownames(sample_data) %in%
                                            common_ids, , drop = FALSE])
    X_matched <- X[rownames(X) %in% common_ids, , drop = FALSE]

    if (!identical(rownames(sdata_matched), rownames(X_matched))) {
      # Attempt to reorder X_matched to match sdata_matched
      X_matched <- X_matched[rownames(sdata_matched), , drop = FALSE]
      # Re-check alignment
      if (!identical(rownames(sdata_matched), rownames(X_matched))) {
        stop("Failed to align rownames of 'X' and 'sample_data'. Please check sample identifiers.")
      }
    }
  } else {
    X_matched <- X
    sdata_matched <- sample_data
  }
  return(list(X = X_matched,
              sample_data = sdata_matched))
}

#' @title Remove samples from a matrix based on IDs or pattern
#'
#' @description
#' Removes specified samples from the input matrix `X` either by providing
#' explicit sample identifiers or a regular expression pattern. Optionally,
#' warns when any provided identifiers are not found.
#'
#' @param X A matrix with samples in rows and features in columns. Row names
#'  must correspond to sample identifiers.
#' @param remove_ids A character vector of sample IDs to remove, or a single
#'  string used as a regular expression pattern. If `NULL`, no samples are removed.
#'
#' @return A matrix with the specified samples removed.
#'
#' @keywords internal
#'
#' @noRd
remove_samples <- function(X, remove_ids = NULL) {
  if (is.null(remove_ids)) {
    X_filtered <- X
  } else {
    if (!is.character(remove_ids)) {
      warning("'remove_ids' will be coerced to character type.\n",
              immediate. = TRUE)
      remove_ids <- as.character(remove_ids)
    }
    if (length(remove_ids) == 1) {
      discard <- grep(remove_ids, rownames(X))
      if (length(discard) == 0) {
        stop("Pattern in 'remove_ids' not found in rownames(X).")
      }
    } else {
      discard <- which(rownames(X) %in% remove_ids)
      unmatched <- setdiff(remove_ids, rownames(X))
      if (length(unmatched) > 0) {
        warning("The following sample IDs were not found in 'X': ",
                paste(unmatched, collapse = ", "))
      }
    }
    X_filtered <- X[-discard, , drop = FALSE]
    if (nrow(X_filtered) == 0) {
      stop("No samples remain after filtering by 'remove_ids'")
    }
  }
  return(X_filtered)
}

#' @title Remove features with low prevalence
#'
#' @description
#' Filters features (columns) from a matrix based on minimum prevalence
#' across samples. Missing values can be defined as either `NA` or `0`.
#'
#' @param X A numeric matrix with samples in rows and features in columns.
#' @param label Character value indicating what constitutes a missing value:
#'  `"NA"` (default) or `"0"`. Only one option should be provided.
#' @param min_prev Numeric value in (0, 1], the minimum proportion of
#'  non-missing samples required to retain a feature.
#' @param verbose Logical. If `TRUE`, prints a message with the number
#'  of features retained.
#'
#' @return A filtered matrix `X` with features removed based on prevalence.
#'
#' @keywords internal
#'
#' @noRd
remove_features <- function(X, label = c("NA", "0"), min_prev,
                            verbose = TRUE) {
  label <- match.arg(label)

  if (!is.numeric(min_prev) || min_prev <= 0 || min_prev > 1) {
    stop("'min_prev' must be a numeric value > 0 and <= 1.")
  }

  if (label == "NA") {
    missing_values <- colSums(is.na(X))
  } else if (label == "0") {
    missing_values <- colSums(X == 0)
  }

  keep_cols <- which(missing_values <= (1 - min_prev) * nrow(X))
  if (length(keep_cols) == 0) {
    stop("No features passed the 'min_prev' theshold. Adjust 'min_prev'.")
  } else {
    X_prev <- X[, keep_cols, drop = FALSE]
    if (verbose) {
      message(ncol(X_prev), " out of ", ncol(X), " features were kept after " ,
              100*min_prev, " % prevalence filter.\n")
    }
  }
  return(X_prev)
}

#' @title Check and format NGS-derived tables
#'
#' @description
#' Internal helper function for verifying consistency between the NGS feature
#' matrix `X`, a taxonomy table, and an optional phylogenetic tree.
#' It ensures the dimensions and names match, shortens feature IDs if necessary,
#' and removes stored DNA sequence columns to avoid issues when constructing a phyloseq object.
#'
#' @param X A numeric matrix with samples in rows and NGS features (e.g., ASVs) in columns.
#' @param taxa_table Optional. A data frame containing taxonomy annotations.
#'  Row names must match `colnames(X)`.
#' @param phylo_tree Optional. A phylogenetic tree.
#'
#' @return A list with:
#' \describe{
#'   \item{X}{The input feature matrix.}
#'   \item{taxa_table}{Formatted taxonomy table with harmonized IDs and cleaned columns.}
#' }
#'
#' @keywords internal
#'
#' @noRd
check_ngs_input <- function(X, taxa_table, phylo_tree, verbose = TRUE) {
  if (is.null(taxa_table) && verbose) {
    message("The phyloseq object will be built without taxonomy table: 'taxa_table' not provided.\n")
  } else {
    # Check that feature ids match
    if (dim(X)[2] != dim(taxa_table)[1]) {
      stop("Number of features in `X` columns and `taxa_table` rows should be equal.")
    } else if (!identical(colnames(X), rownames(taxa_table))) {
      stop("Feature colnames in 'X' and 'taxa_table' rownames should be identical.")
    } else {
      # Check for long ids (e.g. DNA sequences)
      long_ids <- max(nchar(rownames(taxa_table)))
      if (long_ids > 10) {
        warning("Feature names are too long (are they DNA sequences?).\n",
                "Shorter placeholder ids should be used to improve computer memory and speed.\n",
                immediate. = TRUE)
      }
    }
    # Format taxa_table columns
    colnames(taxa_table) <- tolower(colnames(taxa_table))
    sequence_col <- grep("sequence", colnames(taxa_table))
    if (length(sequence_col) != 0) {
      warning("Column named 'sequence' or 'sequences' was detected. ",
              "It will be removed from 'taxa_table' to avoid unwanted phyloseq behaviour.\n",
              immediate. = TRUE)
      taxa_table <- taxa_table[, -sequence_col, drop = FALSE]
    }
    taxa_table <- as.matrix(taxa_table)
  }
  if (is.null(phylo_tree) && verbose) {
    message("The phyloseq object will be built without phylogenetic tree: 'phylo_tree' not provided.\n")
  }
  return(list(X = X,
              taxa_table = taxa_table))
}
