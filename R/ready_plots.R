#' @title Generate plots from a differential analysis (`dana`) object
#'
#' @description
#' This function produces a range of coefficient- and feature-level plots from a `dana` object
#' for a given model term of interest. It supports both main effect and interaction terms,
#' and can visualize significant results from either `fit` or `lrt` P values.
#'
#' @details
#' When `add_interactions = TRUE`, the function shows `fit` coefficients that
#' match significant main and interaction terms.
#'
#' If no significant features are found under the specified `alpha` significance
#' threshold, the function will abort.
#'
#' @param dana_obj A `dana` object returned by `dana()`, containing model results.
#' @param term_name The name of the model term to plot (e.g., `"group"` or `"group:time"`).
#' @param pval_match Regex pattern to match the desired P value column in the results.
#' @param alpha Numeric. Significance threshold to consider features for plotting. Default 0.1.
#' @param add_interactions Logical. Whether to include interaction terms related to `term_name`.
#' @param add_labels Logical. Whether to add custom feature labels in plots.
#'  A "feat_name" or "taxon_name" column must be in the `dana` object.
#'  See `add_taxa()` and `add_feat_name()`.
#' @param plot_coeff Logical. Whether to generate coefficient-level plots.
#'  Will generate volcano, heatmap and dot plots.
#' @param plot_feat Logical. Whether to generate feature-level plots for a specific
#'  variable in `sample_data`.
#' @param plot_ranef Logical. Whether to generate random effect variance plots.
#'  Only for mixed-effects models.
#' @param X_colnames Optional. Character vector specifying which features from `X` to plot.
#'  If `NULL` and `plot_feat = TRUE` (the default), top 10 features based on P value are selected.
#' @param sdata_var Character. A column in `dana_obj$sdata` used for feature-level
#'  plots when `plot_feat = TRUE`.
#' @param group_colours Optional named vector of colours for `sdata_var` groups
#'  to be passed as `values` argument to `ggplot2::scale_fill_manual()`.
#' @param paired_id Optional. Column name in `sdata` specifying sample pairing (e.g., subject_id).
#' @param verbose Logical. Whether to display messages during processing.
#' @param ... Additional `ggplot2::theme()` arguments passed to internal plotting helpers (e.g., font sizes).
#'
#' @return A named list of `ggplot` objects stored in `dana_obj$plots`. These may include:
#' \itemize{
#'   \item \code{coeff_volcano}, \code{coeff_heatmap}, \code{coeff_point}
#'   \item \code{feat_scatter}, \code{feat_boxplot}, \code{feat_violin}, \code{feat_ridge}
#'   \item \code{ranef_all}
#' }
#'
#' @seealso
#' * [dana()] for fitting differential analysis models on omics datasets.
#' * [add_taxa()] and [add_feat_name()] for adding feature labels to dana object.
#' * [ggplot2::ggplot()] and [ggplot2::theme()] to further customise plots.
#'
#' @examples
#' set.seed(123)
#' mock_X <- matrix(rnorm(20 * 5), nrow = 20)
#' colnames(mock_X) <- paste0("feat_", seq_len(5))
#' rownames(mock_X) <- paste0("sample_", seq_len(20))
#'
#' sample_data <- data.frame(
#'   sample_id = rownames(mock_X),
#'   group = factor(rep(c("A", "B"), each = 10)),
#'   time = factor(rep(c("T1", "T2"), times = 10)),
#'   subject_id = factor(rep(seq_len(10), each = 2)),
#'   stringsAsFactors = FALSE
#' )
#' rownames(sample_data) <- sample_data$sample_id
#'
#' fit_df <- data.frame(
#'   feat_id = rep(colnames(mock_X), each = 2),
#'   Coefficient = rep(c("(Intercept)", "groupB"), 5),
#'   Estimate = rnorm(10),
#'   `Pr(>|t|)` = runif(10),
#'   padj = runif(10),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Mock dana object
#' dana_obj <- list(
#'   X = mock_X,
#'   sdata = sample_data,
#'   formula_rhs = ~ group,
#'   fit = fit_df,
#'   lrt = data.frame(),  #' empty but valid
#'   ranef = data.frame() #' empty but valid
#' )
#' class(dana_obj) <- "dana"
#'
#' dana_obj <- dana_obj |>
#' ready_plots(
#'   term_name = "group",
#'   pval_match = "padj",
#'   alpha = 0.5,
#'   add_labels = FALSE,
#'   plot_coeff = TRUE,
#'   plot_feat = TRUE,
#'   plot_ranef = FALSE,
#'   sdata_var = "group",
#'   verbose = FALSE
#' )
#'
#' # Visualize generated plots
#' dana_obj$plots
#'
#' @importFrom rlang .data
#'
#' @export
ready_plots <- function(dana_obj, term_name, pval_match, alpha = 0.1,
                        add_interactions = TRUE, add_labels = TRUE,
                        plot_coeff = TRUE, plot_feat = TRUE, plot_ranef = FALSE,
                        X_colnames = NULL, sdata_var = NULL, group_colours = NULL,
                        paired_id = NULL, verbose = TRUE, ...) {
  # Match or set arguments
  X <- dana_obj$X
  fit_df <- dana_obj$fit
  lrt_df <- dana_obj$lrt
  ranef_df <- dana_obj$ranef
  term_model <- attr(stats::terms(dana_obj$formula_rhs), "term.labels")
  feat_label <- grep("_name", colnames(fit_df), value = TRUE)

  # Check parameters
  if (!plot_coeff && !plot_feat && !plot_ranef) {
    stop("At least one of 'plot_coeff', 'plot_feat', or 'plot_ranef' must be TRUE.")
  }

  if (!term_name %in% term_model) {
    stop("'term_name' must exist in 'formula_rhs' from 'dana' class object.")
  }

  stopifnot("'alpha' must be numeric" = is.numeric(alpha),
            "'alpha' < 0" = alpha > 0,
            "'alpha' > 1" = alpha < 1)

  if (add_labels) {
    if (length(feat_label) == 0) {
      warning("Cannot find 'feat_name' or 'taxon_name' in dana object. ",
              "Plots will show default feature labels.\n",
              immediate. = TRUE)
      feat_label <- "feat_id"
    } else if (length(feat_label) > 1 && all(feat_label %in% c("feat_name", "taxon_name"))) {
      warning("Both 'feat_name' and 'taxon_name' columns found. 'taxon_name' will be used as labels.\n",
              immediate. = TRUE)
      feat_label <- "taxon_name"
    }
  } else {
    feat_label <- "feat_id"
  }

  if (plot_feat) {
    if (is.null(sdata_var) || !is.character(sdata_var)) {
      stop("'sdata_var' must be a character string.")
    } else if (!sdata_var %in% colnames(dana_obj$sdata)) {
      stop("'sdata_var' not found in sample data.")
    }
    if (!is.null(paired_id) && !paired_id %in% colnames(dana_obj$sdata)) {
      stop("'paired_id' variable was not found in sample data.")
    }
    if (!is.null(X_colnames) &&
        (!is.character(X_colnames) ||
         any(!X_colnames %in% colnames(X)))) {
      stop("'X_colnames' must be a character vector of omics data colnames or NULL.")
    }
  }

  if (plot_ranef && nrow(ranef_df) == 0) {
    stop("Cannot find random effects: 'plot-ranef' is only possible for mixed-effects models.")
  }

  # P values and feature filtering
  p_nominal <- grep("Pr", colnames(fit_df), value = TRUE)
  if (length(p_nominal) != 1) {
    stop("Expected exactly one nominal P-value column, found: ", toString(p_nominal))
  }

  padj_colname <- grep(pval_match, colnames(fit_df), value = TRUE)
  if (length(padj_colname) != 1) {
    stop("'pval_match' pattern should match one column of P values. Matched: ",
         paste(padj_colname, collapse = ", "))
  }

  label_term <- unique(grep(term_name, fit_df[["Coefficient"]], value = TRUE))
  has_interactions <- any(grepl(":", label_term))
  all_interactions <- all(grepl(":", label_term))

  keep_term <- grepl(term_name, fit_df[["Coefficient"]])
  if (!add_interactions && all_interactions) {
    stop("Set 'add_interactions = TRUE' to plot any significant interaction: ", label_term)
  } else if (!add_interactions && has_interactions) {
    if (verbose) message("Interaction terms will be ignored. ",
                         "Set 'add_interactions = TRUE' to change this behaviour.\n")
    keep_term <- grepl(term_name, fit_df[["Coefficient"]]) & !grepl(":", fit_df[["Coefficient"]])
  } else if (add_interactions && has_interactions) {
    if (verbose) message("Including significant interactions for term '", term_name,
                         "'. Set 'add_interactions = FALSE' to exclude them.\n")
  }

  fit_term <- fit_df[keep_term, , drop = FALSE]
  is_signif <- fit_term[[padj_colname]] < alpha
  n_signif <- sum(is_signif, na.rm = TRUE)
  if (n_signif == 0 && !plot_ranef) {
      stop("No significant results at selected ", alpha, " significance threshold.\n")
  } else if (n_signif == 0 && plot_ranef) {
    warning("No significant results at selected ", alpha, " significance threshold. ",
            "Only random effects will be plotted.\n")
    plot_coeff <- FALSE
    plot_feat <- FALSE
  } else {
    if (verbose) message(n_signif, " significant results found at selected ", alpha,
                         " significance threshold.\n")
  }

  term_match <- stringr::str_split_1(term_name, pattern = ":")
  term_n <- length(term_match)

  # Generate plots
  plots <- list()

  if (plot_coeff) {
    # Volcano
    plots$coeff_volcano <- plot_volcano(fit_term, padj_colname, alpha, feat_label)

    # Heatmap
    if (n_signif > 50 && verbose) {
      message("Only the first 50 most significant features will be used for the heatmap plot.\n")
    }
    p_top50 <- fit_term |>
      dplyr::filter(.data[[padj_colname]] < alpha) |>
      dplyr::arrange(.data[[padj_colname]],
                     .data[[p_nominal]],
                     dplyr::desc(.data[["Estimate"]])) |>
      dplyr::slice_head(n = 50)

    plots$coeff_heatmap <- plot_heatmap(p_top50, feat_label, ...)

    # Dot plot
    plots$coeff_point <- plot_point(p_top50, padj_colname, feat_label, ...)
  }

  if (plot_feat) {
    # Determine feature IDs to plot
    if (!is.null(X_colnames)) {
      selected_feat_ids <- X_colnames
    } else {
      if (verbose) message("'X_colnames' = NULL: the first 10 most significant features will be plotted.\n")
      selected_feat_ids <- fit_term |>
        dplyr::filter(.data[[padj_colname]] < alpha) |>
        dplyr::arrange(.data[[padj_colname]],
                       .data[[p_nominal]],
                       dplyr::desc(.data[["Estimate"]])) |>
        dplyr::slice_head(n = 10) |>
        dplyr::pull(.data[["feat_id"]])
    }

    # Subset omic data matrix
    missing_feats <- setdiff(selected_feat_ids, colnames(X))
    if (length(missing_feats) > 0) {
      warning("Some selected features not found in X: ", paste(missing_feats, collapse = ", "))
    }
    X_selected <- X[, which(colnames(X) %in% selected_feat_ids), drop = FALSE]

    # Prepare long format
    p_selected <- data.frame(sample_id = rownames(X_selected),
                             X_selected,
                             check.names = FALSE) |>
                  tidyr::pivot_longer(-"sample_id",
                                      names_to = "feat_id",
                                      values_to = "value")
    p_selected <- dplyr::left_join(p_selected,
                                   fit_term[, c("feat_id", feat_label)],
                                   by = "feat_id")

    p_feat <- dplyr::left_join(p_selected,
                               dana_obj$sdata,
                               by = "sample_id")

    if (is.numeric(p_feat[[sdata_var]])) {
      # Scatter plot
      plots$feat_scatter <- plot_feat(p_feat, plot_type = "scatter", sdata_var, paired_id,
                                      feat_label, group_colours = NULL, ...)
    } else {
      # Boxplot
      plots$feat_boxplot <- plot_feat(p_feat, plot_type = "boxplot", sdata_var, paired_id,
                                      feat_label, group_colours, ...)
      # Violin
      plots$feat_violin <- plot_feat(p_feat, plot_type = "violin", sdata_var, paired_id,
                                     feat_label, group_colours, ...)
      # Ridgeline
      plots$feat_ridge <- plot_feat_ridge(p_feat, sdata_var, feat_label, group_colours, ...)
    }
  }

  if (plot_ranef) {
    p_ranef <- ranef_df |> dplyr::rename(value = .data[["vcov"]])

    plots$ranef_all <- ggplot2::ggplot(data = p_ranef,
                                       ggplot2::aes(x = .data[["grp"]],
                                                    y = .data[["value"]],
                                                    fill = .data[["grp"]])) +
      ggplot2::geom_violin(draw_quantiles = TRUE,
                           show.legend = FALSE,
                           linewidth = 0.05,
                           scale = "width") +
      ggplot2::geom_boxplot(width = 0.1,
                            outlier.size = 0.25,
                            linewidth = 0.05) +
      ggplot2::scale_fill_brewer(palette = "BrBG") +
      ggplot2::labs(x = NULL,
                    y = "Variance") +
      ready_theme(...)
  }

  dana_obj$plots <- plots

  return(dana_obj)
}
