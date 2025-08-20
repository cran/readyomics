#' @title Multivariate analysis (PCA, PLS, OPLS)
#'
#' @description
#' Performs PCA, PLS, or OPLS using ropls and generates a formatted scores plot
#' based on the first two components.
#'
#' @details
#' The analysis type depends on the `...` arguments passed to `ropls::opls()`.
#'
#' @param X A numeric matrix or data frame of features (e.g., metabolites, genes),
#'   with samples as rows and features as columns.
#' @param sample_data A `data.frame` containing sample-level data. Row names must match
#'   the sample identifiers in `X` and must be also in a column named `"sample_id"`.
#' @param group_colour Optional. Character colname in `sample_data` used for point color mapping.
#' @param group_shape Optional. Character colname in `sample_data` used for point shape mapping.
#' @param plot_title Optional. Character string specifying the plot title.
#' @param verbose Logical. If `TRUE`, displays progress messages.
#' @param ... Additional arguments passed to `ropls::opls()` (e.g.`predI =`, `orthoI =`).
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{ropls_obj}{The `ropls::opls()` object.}
#'   \item{scores_plot}{A `ggplot2::ggplot()` object showing the scores plot.}
#' }
#'
#' @seealso [ropls::opls()] for details on the `ropls::opls()` output.
#'
#' @examples
#' # PCA
#' set.seed(123)
#' mock_X <- matrix(rnorm(40),
#'                  nrow = 10,
#'                  dimnames = list(paste0("sample", 1:10),
#'                                  paste0("feat", 1:4))
#'                  )
#'
#' sample_data <- data.frame(
#'   sample_id = rownames(mock_X),
#'   group = factor(rep(c("A", "B"), each = 5)),
#'   batch = factor(rep(1:2, times = 5)),
#'   row.names = rownames(mock_X),
#'   stringsAsFactors = FALSE
#' )
#'
#' result <- mva(
#'   X = mock_X,
#'   sample_data = sample_data,
#'   group_colour = "group",
#'   group_shape = "batch",
#'   plot_title = "Test PCA Plot",
#'   predI = 2,  # PCA: set components
#'   verbose = FALSE
#' )
#'
#' # PCA plot
#' result$scores_plot
#'
#' @importFrom rlang .data
#'
#' @export
mva <- function(X, sample_data, group_colour = NULL,
                group_shape = NULL, plot_title = NULL,
                verbose = TRUE, ...) {
  # Match X and sample_data by sample id
  match_ids <- check_sample_ids(X, sample_data, do_match = TRUE,
                                verbose = verbose)
  X_matched <- match_ids$X
  sdata_matched <- match_ids$sample_data

  # Check parameters
  if (!is.null(group_colour) &&
      (!is.character(group_colour) || !group_colour %in% c(colnames(sample_data)))) {
    stop("Cannot find 'group_colour' in 'sample_data' colnames.")
  }

  if (!is.null(group_shape) &&
      (!is.character(group_shape) || !group_shape %in% c(colnames(sample_data)))) {
    stop("Cannot find 'group_shape' in 'sample_data' colnames.")
  }

  # Run PCA/PLS/OPLS
  mva_obj <- ropls::opls(X_matched, ...)

  # Generate formatted scores plot
  scores <- as.data.frame(mva_obj@scoreMN)
  scores <- cbind(sdata_matched, scores)

  PC1_var <- 100 * mva_obj@modelDF$R2X[1]
  PC2_var <- 100 * mva_obj@modelDF$R2X[2]

  # Set ggplot aesthetics
  p <- ggplot2::ggplot(data = scores,
                       ggplot2::aes(x = .data[["p1"]],
                                    y = .data[["p2"]]))
  if (!is.null(group_colour)) p <- p + ggplot2::aes(colour = .data[[group_colour]])
  if (!is.null(group_shape))  p <- p + ggplot2::aes(shape = .data[[group_shape]])

  p <- p +
    ggplot2::geom_point(size = 0.25,
                        alpha = 0.75) +
    ggplot2::geom_hline(yintercept = 0,
                        linetype = 3,
                        linewidth = 0.05) +
    ggplot2::geom_vline(xintercept = 0,
                        linetype = 3,
                        linewidth = 0.05) +
    ggplot2::labs(x = bquote(t[1] ~ " ("*.(PC1_var)*"%)"),
                  y = bquote(t[2] ~ " ("*.(PC2_var)*"%)"),
                  title = plot_title,
                  shape = NULL,
                  colour = NULL) +
    ggplot2::xlim(-1.1 * max(abs(scores[["p1"]])), 1.1 * max(abs(scores[["p1"]]))) +
    ggplot2::ylim(-1.1 * max(abs(scores[["p2"]])), 1.1 * max(abs(scores[["p2"]]))) +
    ggplot2::theme(axis.title.y = ggplot2::element_text(margin = ggplot2::margin(0, 0, 0, 0)),
                   legend.background = ggplot2::element_rect(fill = NA),
                   legend.key = ggplot2::element_blank(),
                   legend.key.size = grid::unit(0.1, "cm"),
                   legend.margin = ggplot2::margin(0, 0, 0, 0),
                   legend.position = "inside",
                   legend.position.inside = c(0.85, 0.12),
                   legend.spacing = grid::unit(0.1, "cm"),
                   legend.text = ggplot2::element_text(vjust = 0.5,
                                                       size = 4,
                                                       margin = ggplot2::margin(0, 0, 1, 0)),
                   legend.title = ggplot2::element_text(margin = ggplot2::margin(0, 0, 0, 0),
                                                        size = 5),
                   line = ggplot2::element_line(linewidth = 0.05),
                   panel.background = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(fill = NA,
                                                        linetype = "solid",
                                                        linewidth = 0.05),
                   plot.title = ggplot2::element_text(margin = ggplot2::margin(0, 0, 0, 0)),
                   text = ggplot2::element_text(size = 6))

  return(list(ropls_obj = mva_obj,
              scores_plot = p))
}
