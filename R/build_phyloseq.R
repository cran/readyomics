#' @title Build phyloseq objects for all taxonomy ranks
#'
#' @description
#' Constructs a list of phyloseq objects from a feature matrix (`X`), sample data,
#' taxonomy and (optionally) phylogenetic tree data.
#'
#' @details
#' Phyloseq objects for higher taxonomic ranks are also generated when `taxa_table` is provided.
#' Higher rank taxa with labels matching "unclass" or "unknown" are excluded
#' after aggregation.
#'
#' If very long strings are detected as feature IDs in `X` matrix or `taxa_table`,
#' (for example when actual DNA sequence is used as ID), it will issue a warning,
#' as this could significantly slow down computation and increase memory usage.
#'
#' @param X A numeric matrix of NGS features (e.g., ASVs), with samples in rows
#'   and features in columns (recommended) or vice versa.
#' @param sample_data A `data.frame` containing sample data.
#'   Row names must match sample identifiers in `X`.
#' @param taxa_table (Optional) A taxonomy table with row names corresponding to
#'   feature names in `X`, and taxonomic ranks as columns.
#' @param phylo_tree (Optional) A phylogenetic tree.
#' @param taxa_in_rows Logical. If `TRUE`, `X` is assumed to have taxa as rows
#'  and samples as columns.
#' @param verbose Logical. If `TRUE`, diagnostic messages will be printed.
#'
#' @return A named `list` of `phyloseq` objects and related output:
#' \describe{
#'   \item{asv}{Phyloseq object with the raw feature counts (usually ASVs).}
#'   \item{<tax_rank>}{Phyloseq objects of higher taxonomy ranks from `taxa_table`.}
#' }
#'
#' @seealso [phyloseq::phyloseq()] for further details on phyloseq objects.
#'
#' @examples
#' mock_X <- matrix(c(10, 0, 5, 3, 1, 7),
#'                  nrow = 2, byrow = TRUE,
#'                  dimnames = list(c("sample1", "sample2"),
#'                                  c("ASV1", "ASV2", "ASV3"))
#'                  )
#'
#' mock_sample_data <- data.frame(sample_id = c("sample1", "sample2"),
#'                                group = c("A", "B"),
#'                                row.names = c("sample1", "sample2")
#'                                )
#'
#' mock_taxa_table <- data.frame(Domain = c("Bacteria", "Bacteria", "Bacteria"),
#'                               Genus = c("GenusA", "GenusB", "Unknown"),
#'                               row.names = c("ASV1", "ASV2", "ASV3")
#'                               )
#'
#' phyloseq_ready <- build_phyloseq(X = mock_X,
#'                                  sample_data = mock_sample_data,
#'                                  taxa_table = mock_taxa_table,
#'                                  taxa_in_rows = FALSE,
#'                                  verbose = FALSE)
#'
#' @export
build_phyloseq <- function(X, sample_data, taxa_table = NULL,
                           phylo_tree = NULL, taxa_in_rows,
                           verbose = TRUE) {
  # Check X
  if (taxa_in_rows && verbose) {
    message("'taxa_in_rows = TRUE': 'X' will be transposed, with samples in rows.\n")
    X <- t(X) # Features in columns and samples in rows
    taxa_in_rows <- FALSE
  }

  X <- check_omics_matrix(X, platform = "ngs")

  # Check parameters
  check_ngs <- check_ngs_input(X, taxa_table, phylo_tree, verbose = verbose)
  X <- check_ngs$X
  taxa_table <- check_ngs$taxa_table

  # Match X and sample_data by sample id
  match_ids <- check_sample_ids(X, sample_data, do_match = TRUE,
                                verbose = verbose)
  X_matched <- match_ids$X
  sdata_matched <- match_ids$sample_data

  empty_rows <- which(rowSums(X_matched) == 0)
  if (length(empty_rows) > 0 ) {
    X_matched <- X_matched[-empty_rows, , drop = FALSE]
    sdata_matched <- sdata_matched[-empty_rows, , drop = FALSE] |> droplevels()
  }
  stopifnot(nrow(X_matched) > 0)

  empty_cols <- which(colSums(X_matched) == 0)
  if (length(empty_cols) > 0 ) {
    X_matched <- X_matched[, -empty_cols, drop = FALSE]
    sdata_matched <- sdata_matched[, -empty_cols, drop = FALSE] |> droplevels()
  }
  stopifnot(ncol(X_matched) > 0)

  # Build phyloseq objects
  phyloseq_obj <- list()

  phyloseq_obj$asv <- phyloseq::phyloseq(
    phyloseq::otu_table(X_matched,
                        taxa_are_rows = taxa_in_rows),
    phyloseq::tax_table(taxa_table, errorIfNULL = FALSE),
    phyloseq::sample_data(sdata_matched),
    phylo_tree)

  if (!is.null(taxa_table)) {
    taxa_levels <- colnames(taxa_table)
    highest_rank <- taxa_levels[1] # Assuming descending taxonomic hierarchy
    highest_n <- length(unique(as.character(taxa_table[, highest_rank])))
    if (highest_n == 1) taxa_levels <- taxa_levels[-1] # Single Domain tables
    for (i in taxa_levels) {
      phyloseq_temp <- phyloseq::tax_glom(phyloseq_obj$asv, taxrank = i)
      unclass_taxa <- grep("unclass|Unknown", as.data.frame(phyloseq_temp@tax_table)[[i]])
      if (length(unclass_taxa) > 0) {
        assigned_taxa <- phyloseq::taxa_names(phyloseq_temp)[-unclass_taxa]
        phyloseq_temp <- phyloseq::prune_taxa(assigned_taxa, phyloseq_temp)
      }
      phyloseq_obj[[i]] <- phyloseq_temp
    }
  }

  return(phyloseq_obj)
}
