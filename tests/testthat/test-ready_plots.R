test_that("ready_plots generates plots from a minimal dana object (fit analysis)", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("stringr")

  set.seed(123)
  mock_X <- matrix(rnorm(20 * 5), nrow = 20)
  colnames(mock_X) <- paste0("feat_", seq_len(5))
  rownames(mock_X) <- paste0("sample_", seq_len(20))

  sample_data <- data.frame(
    sample_id = rownames(mock_X),
    group = rep(c("A", "B"), each = 10),
    time = rep(c("T1", "T2"), times = 10),
    subject_id = rep(seq_len(10), each = 2),
    stringsAsFactors = FALSE
  )
  rownames(sample_data) <- sample_data$sample_id

  fit_df <- data.frame(
    feat_id = rep(colnames(mock_X), each = 2),
    Coefficient = rep(c("(Intercept)", "groupB"), 5),
    Estimate = rnorm(10),
    `Pr(>|t|)` = runif(10),
    padj = runif(10),
    stringsAsFactors = FALSE
  )

  dana_obj <- list(
    X = mock_X,
    sdata = sample_data,
    formula_rhs = ~ group,
    fit = fit_df,
    lrt = data.frame(),  # empty but valid
    ranef = data.frame() # empty but valid
  )
  class(dana_obj) <- "dana"

  result <- ready_plots(
    dana_obj = dana_obj,
    term_name = "group",
    pval_match = "padj",
    alpha = 0.5,
    add_labels = FALSE,
    plot_coeff = TRUE,
    plot_feat = TRUE,
    plot_ranef = FALSE,
    sdata_var = "group",
    verbose = FALSE
  )

  expect_type(result, "list")
  expect_true("coeff_volcano" %in% names(result$plots))
  expect_s3_class(result$plots$coeff_volcano, "ggplot")
  expect_true("coeff_heatmap" %in% names(result$plots))
  expect_s3_class(result$plots$feat_boxplot, "ggplot")
})

test_that("ready_plots throws error if no significant features are found", {
  dana_obj <- list(
    X = matrix(rnorm(20), ncol = 2,
               dimnames = list(paste0("s", 1:10), c("feat1", "feat2"))),
    sdata = data.frame(
      sample_id = paste0("s", 1:10),
      group = rep(c("A", "B"), each = 5),
      stringsAsFactors = FALSE,
      row.names = paste0("s", 1:10)
    ),
    formula_rhs = ~ group,
    fit = data.frame(
      feat_id = rep(c("feat1", "feat2"), each = 2),
      Coefficient = rep(c("(Intercept)", "groupB"), 2),
      Estimate = rnorm(4),
      padj = rep(1, 4),
      `Pr(>|t|)` = runif(4),
      stringsAsFactors = FALSE
    ),
    lrt = data.frame(),
    ranef = data.frame()
  )
  class(dana_obj) <- "dana"

  expect_error(
    ready_plots(
      dana_obj,
      term_name = "group",
      pval_match = "padj",
      alpha = 0.05,
      add_labels = FALSE,
      sdata_var = "group",
      verbose = FALSE
    ),
    "No significant results"
  )
})

test_that("ready_plots throws error with missing sdata_var", {
  dana_obj <- list(
    X = matrix(rnorm(20), ncol = 2,
               dimnames = list(paste0("s", 1:10), c("f1", "f2"))),
    sdata = data.frame(
      sample_id = paste0("s", 1:10),
      group = rep(c("A", "B"), each = 5),
      stringsAsFactors = FALSE,
      row.names = paste0("s", 1:10)
    ),
    formula_rhs = ~ group,
    fit = data.frame(
      feat_id = rep("f1", 2),
      Coefficient = c("(Intercept)", "groupB"),
      Estimate = rnorm(2),
      padj = c(0.05, 0.02),
      `Pr(>|t|)` = c(0.04, 0.01),
      stringsAsFactors = FALSE
    ),
    lrt = data.frame(),
    ranef = data.frame()
  )
  class(dana_obj) <- "dana"

  expect_error(
    ready_plots(
      dana_obj,
      term_name = "group",
      pval_match = "padj",
      alpha = 0.1,
      add_labels = FALSE,
      plot_feat = TRUE,
      sdata_var = "missing_column",
      verbose = FALSE
    ),
    "'sdata_var' not found in sample data"
  )
})

test_that("ready_plots supports paired_id if present", {
  skip_if_not_installed("ggplot2")

  mock_dana <- list(
    X = matrix(rnorm(20 * 5), ncol = 5,
               dimnames = list(paste0("sample", 1:20), paste0("feat", 1:5))),
    sdata = data.frame(
      sample_id = paste0("sample", 1:20),
      group = rep(c("A", "B"), each = 10),
      subject_id = rep(1:10, each = 2),
      stringsAsFactors = FALSE,
      row.names = paste0("sample", 1:20)
    ),
    formula_rhs = ~ group,
    fit = data.frame(
      feat_id = rep(paste0("feat", 1:5), each = 2),
      Coefficient = rep(c("(Intercept)", "groupB"), 5),
      Estimate = rnorm(10),
      padj = runif(10),
      `Pr(>|t|)` = runif(10),
      stringsAsFactors = FALSE
    ),
    lrt = data.frame(),
    ranef = data.frame()
  )
  class(mock_dana) <- "dana"

  result <- ready_plots(
    dana_obj = mock_dana,
    term_name = "group",
    pval_match = "padj",
    alpha = 0.5,
    add_labels = FALSE,
    plot_feat = TRUE,
    plot_ranef = FALSE,
    sdata_var = "group",
    paired_id = "subject_id",
    verbose = FALSE
  )

  expect_type(result, "list")
  expect_true("feat_boxplot" %in% names(result$plots))
  expect_s3_class(result$plots$feat_boxplot, "ggplot")
})

test_that("ready_plots falls back to 'feat_id' when no label columns are found", {
  dana_obj <- list(
    X = matrix(rnorm(20 * 4), ncol = 4,
               dimnames = list(paste0("sample", 1:20), paste0("f", 1:4))),
    sdata = data.frame(
      sample_id = paste0("sample", 1:20),
      group = rep(c("A", "B"), each = 10),
      row.names = paste0("sample", 1:20)
    ),
    formula_rhs = ~ group,
    fit = data.frame(
      feat_id = rep(paste0("f", 1:4), each = 2),
      Coefficient = rep(c("(Intercept)", "groupB"), 4),
      Estimate = rnorm(8),
      padj = runif(8),
      `Pr(>|t|)` = runif(8),
      stringsAsFactors = FALSE
    ),
    lrt = data.frame(),
    ranef = data.frame()
  )
  class(dana_obj) <- "dana"

  expect_warning({
    result <- ready_plots(
      dana_obj = dana_obj,
      term_name = "group",
      pval_match = "padj",
      alpha = 0.5,
      plot_feat = TRUE,
      sdata_var = "group",
      verbose = FALSE
    )
    "Default labels will be used"
  })

  expect_true("feat_boxplot" %in% names(result$plots))
  expect_s3_class(result$plots$feat_boxplot, "ggplot")
})

test_that("ready_plots respects custom X_colnames input", {
  set.seed(1)
  dana_obj <- list(
    X = matrix(rnorm(100), nrow = 10,
               dimnames = list(paste0("s", 1:10), paste0("f", 1:10))),
    sdata = data.frame(
      sample_id = paste0("s", 1:10),
      group = rep(c("A", "B"), each = 5),
      stringsAsFactors = FALSE,
      row.names = paste0("s", 1:10)
    ),
    formula_rhs = ~ group,
    fit = data.frame(
      feat_id = rep(paste0("f", 1:10), each = 2),
      Coefficient = rep(c("(Intercept)", "groupB"), 10),
      Estimate = rnorm(20),
      padj = rep(0.01, 20),
      `Pr(>|t|)` = runif(20),
      stringsAsFactors = FALSE
    ),
    lrt = data.frame(),
    ranef = data.frame()
  )
  class(dana_obj) <- "dana"

  result <- ready_plots(
    dana_obj = dana_obj,
    term_name = "group",
    pval_match = "padj",
    alpha = 0.05,
    add_labels = FALSE,
    plot_feat = TRUE,
    sdata_var = "group",
    X_colnames = c("f1", "f3", "f5"),
    verbose = FALSE
  )

  expect_true("feat_boxplot" %in% names(result$plots))
  expect_s3_class(result$plots$feat_boxplot, "ggplot")
})

test_that("ready_plots throws error if term_name is not in model", {
  dana_obj <- list(
    X = matrix(rnorm(100), nrow = 10,
               dimnames = list(paste0("s", 1:10), paste0("f", 1:10))),
    sdata = data.frame(
      sample_id = paste0("s", 1:10),
      group = rep(c("A", "B"), each = 5),
      row.names = paste0("s", 1:10)
    ),
    formula_rhs = ~ group,
    fit = data.frame(
      feat_id = rep(paste0("f", 1:10), each = 2),
      Coefficient = rep(c("(Intercept)", "groupB"), 10),
      Estimate = rnorm(20),
      padj = runif(20),
      `Pr(>|t|)` = runif(20),
      stringsAsFactors = FALSE
    ),
    lrt = data.frame(),
    ranef = data.frame()
  )
  class(dana_obj) <- "dana"

  expect_error(
    ready_plots(
      dana_obj = dana_obj,
      term_name = "not_a_term",
      pval_match = "padj",
      alpha = 0.1,
      plot_feat = TRUE,
      sdata_var = "group",
      verbose = FALSE
    ),
    "'term_name' must exist in 'formula_rhs'"
  )
})
