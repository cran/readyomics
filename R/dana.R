#' @title Differential analysis (dana)
#'
#' @description
#' Feature-wise `stats::lm()` or `lme4::lmer()` models of an omics data matrix.
#' Supports likelihood ratio tests (LRT) and parallel computation.
#'
#' @details
#' Models are fit independently for each feature using `stats::lm()` or `lmerTest::lmer()`,
#' depending on whether `dana()` detects random effects in `formula_rhs`.
#' Feature-wise models can be evaluated in parallel using `future::plan()`,
#' with optional progress updates via `progressr::with_progress()`.
#'
#' @param X A numeric matrix with samples in rows and features in columns. Sample
#'  IDs in row names must match the format from `sample_id` column in `sample_data`.
#' @param sample_data A data frame containing sample-level data.
#'   Must have a `sample_id` column matching row names in `X` and `sample_data`.
#' @param formula_rhs A one-sided formula (e.g., `~ group + (1|subject)`).
#'   Must not contain a response variable.
#' @param term_LRT Optional. Character vector of formula terms to test via LRT.
#'  Random effects must be written without parentheses (e.g., `"1 | group"`).
#' @param model_control Optional. List of control arguments passed to the model.
#' @param platform Character string indicating the omics platform (e.g., `"ms"`, `"nmr"`, `"ngs"`).
#' @param assay Optional. Character string indicating the name of the platform assay (e.g., `"lipidomics"`).
#' @param verbose Logical. If TRUE, prints progress messages.
#'
#' @return An object of class `"dana"`:
#' \describe{
#'   \item{X}{Matched data matrix.}
#'   \item{sdata}{Matched sample data.}
#'   \item{fit}{Data frame of model coefficients and confidence intervals per feature.}
#'   \item{lrt}{Likelihood ratio test results (if `term_LRT` is specified).}
#'   \item{ranef}{Random effects variance components (if using mixed models).}
#'   \item{errors}{A data frame logging any model fitting errors per feature.}
#' }
#'
#' @seealso [stats::lm()], [lme4::lmer()], [lmerTest::lmer()] parameters.
#'
#' @examples
#' mock_X <- matrix(
#'   rnorm(50 * 10) +
#'     rep(c(rep(0, 25), rep(2, 25)), each = 10) * rep(1:10 %in% 1:3, each = 50),
#'   nrow = 50
#' )
#'
#' rownames(mock_X) <- paste0("sample", 1:50)
#' colnames(mock_X) <- paste0("feat", 1:10)
#'
#' sample_data <- data.frame(
#'   sample_id = rownames(mock_X),
#'   group = factor(rep(c("A", "B"), each = 25)),
#'   subject = factor(rep(1:25, each = 2)),
#'   row.names = rownames(mock_X)
#' )
#'
#' # Example with parallel computation setup (not run)
#' # future::plan(multisession)
#' # progressr::handlers(global = TRUE)
#' # progressr::with_progress({
#'   result <- dana(X = mock_X,
#'                  sample_data = sample_data,
#'                  formula_rhs = ~ group + (1 | subject),
#'                  term_LRT = c("group", "1 | subject"), # Multiple terms allowed
#'                  platform = "ms",
#'                  assay = "lipidomics",
#'                  verbose = FALSE
#'                  )
#' # })
#'
#' # Modify `dana` object at once with pipes (not run)
#' # dana_obj <- dana_obj |> adjust_pval() |> add_feat_name() |> ready_plots()
#'
#' @export
dana <- function(X, sample_data, formula_rhs, term_LRT = NULL,
                 model_control = list(), platform = c("ms", "nmr", "ngs"),
                 assay = NULL, verbose = TRUE) {
  # Match or set arguments
  platform <- match.arg(platform)

  current_plan <- future::plan()
  if (methods::is(current_plan, "sequential") && ncol(X) > 100) {
    warning("Set 'future::plan' parallel strategy to speed up computation.\n",
            immediate. = TRUE)
  }

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

  if (is.null(term_LRT) && nrow(X) < 50) {
    warning("Likelihood ratio test is recommened for small sample size data (n < 50).\n",
            immediate. = TRUE)
  }

  if (!is.null(term_LRT)) {
    if (!is.character(term_LRT)) {
    stop("'term_LRT' must be a character vector of formula term(s) to perform LRT or NULL.")
    } else if (!all(term_LRT %in% attr(stats::terms(formula_rhs), "term.labels"))) {
      stop("'term_LRT' must be part of the terms included in 'formula_rhs'.")
    }
  }

  if (!is.null(assay) && !is.character(assay)) {
    stop("'assay' must be a character string indicating the name of the omics assay or NULL.")
  }

  # Match X and sample_data by sample id
  match_ids <- check_sample_ids(X, sample_data, do_match = TRUE,
                                verbose = verbose)
  X_matched <- match_ids$X
  sdata_matched <- match_ids$sample_data
  feat_num <- ncol(X_matched)

  # Determine model type
  rand_eff <- lme4::findbars(formula_rhs)
  has_ranef <- length(rand_eff) > 0
  fit_model <- if (has_ranef) lmerTest::lmer else stats::lm
  if (verbose) message("Using model: ", if (has_ranef) "lmerTest::lmer" else "stats::lm.\n")

  # Fit feature-wise models (parallel evaluation recommended)
  if (feat_num < 100) {
    steps <- feat_num
    report <- function(i) TRUE
  } else {
    steps <- floor(feat_num / 10) + 1
    report <- function(i) i %% 10 == 0 || i == feat_num
  }
  p <- progressr::progressor(steps = steps)
  fit_result <- future.apply::future_lapply(seq_len(feat_num), function(i) {
    if (report(i)) p(sprintf("Processed %d of %d", i, feat_num))
    tryCatch({
      data_temp <- cbind(sdata_matched, X_matched[, i, drop = FALSE])
      y_temp <- colnames(X_matched)[i]
      formula_temp <- stats::update(formula_rhs, paste(y_temp, "~ ."))
      model_args_temp <- utils::modifyList(model_control,
                                           list(formula = formula_temp,
                                           data = data_temp))
      model_temp <- do.call(fit_model, model_args_temp)
      coeff_temp <- summary(model_temp)$coefficients
      confint_temp <- if (inherits(model_temp, "lm")) {
                        stats::confint(model_temp)
                      } else {
                        lme4::confint.merMod(model_temp, parm = "beta_")
                      }
      fit_df <- data.frame(Coefficient = rownames(coeff_temp),
                           coeff_temp,
                           confint_temp,
                           feat_id = y_temp,
                           check.names = FALSE,
                           row.names = NULL)

      ranef_varcor <- if (inherits(model_temp, "lm")) {
                          NULL
                        } else {
                          df <- as.data.frame(lme4::VarCorr(model_temp))
                          df$feat_id <- y_temp
                          df
                        }

      lrt_result <- if (is.null(term_LRT)) {
                        NULL
                      } else {
                        lrt_temp <- lapply(term_LRT, function(z) {
                          if (grepl("\\|", z)) { # Random effect
                            term_temp <- paste0("(", z, ")")
                          } else {
                            term_temp <- z
                          }
                          formula_null <- stats::update(formula_temp, paste(". ~ . -", term_temp))
                          model_args_null <- utils::modifyList(model_args_temp,
                                                               list(formula = formula_null))
                          # Test for random effects to choose null model
                          rand_eff_null <- lme4::findbars(formula_null)
                          has_ranef_null <- length(rand_eff_null) > 0
                          fit_model_null <- if (has_ranef_null) lmerTest::lmer else stats::lm
                          model_null <- do.call(fit_model_null, model_args_null)
                          lrt_stats <- if (inherits(model_null, "lm") && inherits(model_temp, "lmerModLmerTest")) {
                                          stats::anova(model_temp, model_null) # Changed order for lmer compatibility
                                        } else {
                                          stats::anova(model_null, model_temp)
                                        }
                          lrt_stats$term <- z
                          lrt_stats
                        })
                        df <- do.call(rbind, lrt_temp)
                        df$feat_id <- y_temp
                        df
                      }
       list(fit = fit_df,
            ranef = ranef_varcor,
            lrt = lrt_result,
            error = NULL)
    }, error = function(e) {
        warning(paste("Model failed for feature", i, ":", e$message))
        list(fit = NULL,
             ranef = NULL,
             lrt = NULL,
             error = data.frame(feat_id = colnames(X_matched)[i], error = e$message))
    })
  }, future.seed = TRUE)

  # Combine output
  fit_all <- data.table::rbindlist(lapply(fit_result, `[[`, "fit"))
  ranef_all <- data.table::rbindlist(lapply(fit_result, `[[`, "ranef"))
  lrt_all <- data.table::rbindlist(lapply(fit_result, `[[`, "lrt"))
  error_log <- data.table::rbindlist(lapply(fit_result, `[[`, "error"))

  # Add platform and assay info if provided
  fit_all$platform <- platform
  if (nrow(ranef_all) != 0) ranef_all$platform <- platform
  if (nrow(lrt_all) != 0) lrt_all$platform <- platform

  if (!is.null(assay)) {
    fit_all$assay <- assay
    if (nrow(ranef_all) != 0) ranef_all$assay <- assay
    if (nrow(lrt_all) != 0) lrt_all$assay <- assay
  }

  return(structure(list(X = X_matched,
                        sdata = sdata_matched,
                        formula_rhs = formula_rhs,
                        fit = fit_all,
                        lrt = lrt_all,
                        ranef = ranef_all,
                        errors = error_log),
                        class = "dana"))
}
