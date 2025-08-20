#' @title PERMANOVA with flexible permutation control
#'
#' @description
#' Performs PERMANOVA (Permutational Multivariate Analysis of Variance).
#' Supports both joint-term (default `vegan::adonis2()`) and single-term testing
#' when `independent = TRUE`.
#' Several distance methods, and fine-grained permutation control.
#'
#' @details
#' - Supports both `stats::dist()` and `vegan::vegdist()` for distance matrix computation.
#' - Distance method must be specified in `dist_control$method`.
#' - Permutation design is controlled via the permute package using `permute::shuffleSet()`.
#' - If `seed` is supplied, the same permutations will be used across runs for reproducibility.
#'
#' @param X A processed matrix or data frame of features (samples in rows, features in columns).
#' @param sample_data A `data.frame` containing sample-level data.
#'   Row names must match those in `X`.
#' @param formula_rhs A one-sided formula (e.g., `~ group + age`).
#' @param dist_control A named list of arguments to control distance calculation.
#'   Must contain at least `method`. Defaults to `"Euclidean"` via `stats::dist()`.
#' @param perm_control A named list specifying `permute::shuffleSet()` parameters.
#'   By default, `joint_terms` parameters will be used, with same `vegan::adonis2()`
#'   defaults, unless variable-specific permutation settings are added as named list
#'   elements (e.g. perm_control = list(joint_terms = , age = , sex = )).
#' @param independent Logical. If `TRUE`, a PERMANOVA test for each variable in `formula_rhs` is performed.
#' @param platform A string specifying the omics platform (`"ms"`, `"nmr"`, `"ngs"`). Used for annotation.
#' @param assay Optional. Character string giving the assay name for annotation (e.g., `"lipidomics"`).
#' @param verbose Logical. If `TRUE`, prints diagnostic messages.
#' @param seed Optional integer. If provided, sets the random seed for reproducible permutation results.
#' @param ... Additional arguments passed to `vegan::adonis2()`.
#'
#' @return A named `list` with three elements:
#' \describe{
#'   \item{X_dist}{A `dist` object.}
#'   \item{perm_matrix_joint}{A `matrix` from `permute::shuffleSet()` joint_terms control.}
#'   \item{permanova_joint}{A `data.frame` of PERMANOVA results using the full model.}
#'   \item{permanova_indep}{A `data.frame` of a PERMANOVA results for each predictor,
#'                          or `NULL` if `independent = FALSE`.}
#' }
#'
#' @seealso
#' * [stats::dist()] and [vegan::vegdist()] for information on available distances.
#' * [vegan::adonis2()] and [permute::shuffleSet()] for control options and details.
#' * [process_ngs()] to pre-process and normalize an `X` NGS dataset.
#' * [process_ms()] to pre-process and normalize an `X` MS dataset.
#'
#' @examples
#' # Mock data
#' X <- matrix(rnorm(40), nrow = 10,
#'             dimnames = list(paste0("sample", 1:10),
#'                             paste0("feat", 1:4)))
#' sample_data <- data.frame(
#'   sample_id = rownames(X),
#'   group = factor(rep(c("A", "B"), each = 5)),
#'   age = rep(20:29, length.out = 10),
#'   row.names = rownames(X),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Simple control structures
#' dist_control <- list(method = "euclidean")
#' perm_control <- list(
#'   joint_terms = list(control = permute::how(blocks = NULL, nperm = 9)),
#'   group = list(control = permute::how(blocks = NULL, nperm = 9)),
#'   age = list(control = permute::how(blocks = NULL, nperm = 9))
#' )
#'
#' result <- permanova(
#'   X = X,
#'   sample_data = sample_data,
#'   formula_rhs = ~ group + age,
#'   dist_control = dist_control,
#'   perm_control = perm_control,
#'   independent = TRUE,
#'   platform = "ms",
#'   assay = "lipidomics",
#'   seed = 42,
#'   verbose = FALSE
#' )
#'
#' @export
permanova <- function(X, sample_data, formula_rhs,
                      dist_control = list(method = "euclidean",
                                          diag = FALSE,
                                          upper = FALSE),
                      perm_control = list(joint_terms = list(control = permute::how(blocks = NULL, nperm = 999))),
                      independent = TRUE, platform = c("ms", "nmr", "ngs"),
                      assay = NULL, seed = NULL, verbose = TRUE, ...) {
  # Match or set arguments
  platform <- match.arg(platform)

  # Check parameters
  if (!inherits(formula_rhs, "formula") ||
      length(formula_rhs) != 2 ||
      attr(stats::terms(formula_rhs), "response") != 0) {
    stop("'formula_rhs' must be a one-sided formula: ~ x1 + x2 + x3")
  }

  var_colnames <- all.vars(formula_rhs)
  missing_vars <- setdiff(var_colnames, colnames(sample_data))
  if (length(missing_vars) > 0) {
    stop("Can't find specified 'formula_rhs' terms in 'sample_data' colnames: ",
         paste(missing_vars, collapse = ", "))
  }

  if (!"joint_terms" %in% names(perm_control)) {
    stop("'perm_control' must contain a 'joint_terms' list with permutation control parameters.")
  }

  missing_perm <- setdiff(var_colnames, names(perm_control))
  if (independent && length(missing_perm) > 0) {
    warning("'independent = TRUE' but no variable-specific permutation controls found. ",
            "Default 'perm_control$joint_terms' will be used.\n",
            immediate. = TRUE)
  }

  if (!is.null(assay) && !is.character(assay)) {
    stop("'assay' must be a character string or NULL.")
  }

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1) {
      stop("'seed' must be a single numeric value.")
    }
    set.seed(seed)
  }

  # Match X and sample_data by sample id
  match_ids <- check_sample_ids(X, sample_data, do_match = TRUE,
                                verbose = verbose)
  X_matched <- match_ids$X
  sdata_matched <- match_ids$sample_data

  # When na.action = na.omit: warn if missing data is removed
  dots <- list(...)
  na.action <- dots$na.action

  if (identical(na.action, stats::na.omit)) {
    sdata_noNA <- stats::model.frame(stats::terms(formula_rhs),
                                     data = sdata_matched,
                                     na.action = na.action)
    sdata_noNA_id <- rownames(sdata_noNA)
    if (nrow(sdata_noNA) < nrow(sdata_matched)) {
      warning("Missing data found in formula_rhs terms: affected rows will be removed.\n",
              immediate. = TRUE)
      match_ids <- check_sample_ids(X, sample_data[which(sample_data$sample_id %in% sdata_noNA_id),],
                                    do_match = TRUE, verbose = verbose)
      X_matched <- match_ids$X
      sdata_matched <- match_ids$sample_data
    }
  }

  # Build distance object
  dist_stats <- c("binary", "canberra", "euclidean", # stats package distances
                  "manhattan", "maximum", "minkowski")
  dist_args <- utils::modifyList(dist_control, list(x = X_matched))

  if (dist_args$method %in% dist_stats) {
    X_dist <- do.call(stats::dist, dist_args)
  } else {
    euclidean_silent <- c("aitchison","chisq","mahalanobis","robust.aitchison")
    if (dist_args$method %in% euclidean_silent && verbose) {
      message(euclidean_silent, "computes Euclidean distance.\n",
              "See 'vegan::vegdist()' documentation and source code for more information.\n")
    }
    X_dist <- do.call(vegan::vegdist, dist_args)
  }

  # Run PERMANOVA
  formula_response <- stats::update(formula_rhs, X_dist ~ .)
  perm_args <- utils::modifyList(perm_control$joint_terms,
                                 list(n = nrow(X_matched)))
  if (!is.null(seed)) set.seed(seed)
  perm_matrix <- do.call(permute::shuffleSet, perm_args)
  permanova_joint <- vegan::adonis2(formula_response, data = sdata_matched,
                                    permutations = perm_matrix, ...)
  permanova_joint$platform <- platform
  permanova_joint$assay <- if (!is.null(assay)) assay else NA

  if (independent) {
    permanova_temp <- list()
    for (i in var_colnames) {
      formula_temp <- stats::as.formula(paste("X_dist ~", i))
      if (i %in% names(perm_control)) {
        # Permutation settings specific for each specified variable
        perm_args_temp <- utils::modifyList(perm_control[[i]],
                                            list(n = nrow(X_matched)))
        if (!is.null(seed)) set.seed(seed)
        perm_matrix_temp <- do.call(permute::shuffleSet, perm_args_temp)
        permanova_temp[[i]] <- vegan::adonis2(formula_temp, data = sdata_matched,
                                              permutations = perm_matrix_temp,
                                              ...)
      } else {
        # Default "joint_terms" parameters will be used otherwise
        if (verbose) {
          message("Permutation parameters for ", i, " have not been specified. ",
                  "Default 'perm_control$joint_terms' will be used.\n")
        }
        permanova_temp[[i]] <- vegan::adonis2(formula_temp, data = sdata_matched,
                                              permutations = perm_matrix,
                                              ...)
      }
      permanova_temp[[i]]$var_name <- i
    }
    permanova_indep <- do.call(rbind, permanova_temp)
    permanova_indep$platform <- platform
    permanova_indep$assay <- if (!is.null(assay)) assay else NA
  } else {
    permanova_indep <- NULL
  }

  return(list(X_dist = X_dist,
              perm_matrix_joint = perm_matrix,
              permanova_joint = permanova_joint,
              permanova_indep = permanova_indep))
}
