test_that("adjust_pval adds BH and storey-adjusted columns to dana_obj$fit", {
  skip_if_not_installed("qvalue")

  set.seed(123)
  feat_ids <- paste0("feat_", 1:4)
  coef_terms <- c("(Intercept)", "groupB")

  dana_obj <- list(
    X = matrix(rnorm(40), nrow = 10, dimnames = list(paste0("s", 1:10), feat_ids)),
    sdata = data.frame(sample_id = paste0("s", 1:10), group = rep(c("A", "B"), each = 5)),
    formula_rhs = ~ group,
    fit = data.frame(
      feat_id = rep(feat_ids, each = 2),
      Coefficient = rep(coef_terms, times = 4),
      Estimate = rnorm(8),
      `Pr(>|t|)` = runif(8),
      stringsAsFactors = FALSE
    ),
    lrt = data.frame(),
    ranef = data.frame()
  )
  class(dana_obj) <- "dana"

  result <- adjust_pval(dana_obj, padj_method = c("BH", "storey"), padj_by = "terms", verbose = FALSE)

  expect_true("padj_BH_groupB" %in% colnames(result$fit))
  expect_true("padj_storey_groupB" %in% colnames(result$fit))
  expect_true(all(is.numeric(result$fit$padj_BH_groupB), na.rm = TRUE))
})

test_that("adjust_pval supports global (all-term) adjustment", {
  set.seed(123)
  feat_ids <- paste0("feat_", 1:4)
  coef_terms <- c("(Intercept)", "groupB")

  dana_obj <- list(
    X = matrix(rnorm(40), nrow = 10, dimnames = list(paste0("s", 1:10), feat_ids)),
    sdata = data.frame(sample_id = paste0("s", 1:10), group = rep(c("A", "B"), each = 5)),
    formula_rhs = ~ group,
    fit = data.frame(
      feat_id = rep(feat_ids, each = 2),
      Coefficient = rep(coef_terms, times = 4),
      Estimate = rnorm(8),
      `Pr(>|t|)` = runif(8),
      stringsAsFactors = FALSE
    ),
    lrt = data.frame(),
    ranef = data.frame()
  )
  class(dana_obj) <- "dana"

  result <- adjust_pval(dana_obj, padj_method = "BH", padj_by = "all", verbose = FALSE)

  expect_true("padj_BH" %in% colnames(result$fit))
  expect_true(all(is.numeric(result$fit$padj_BH), na.rm = TRUE))
})

test_that("adjust_pval falls back to BH if storey fails", {
  skip_if_not_installed("qvalue")
  set.seed(123)
  feat_ids <- paste0("feat_", 1:4)
  coef_terms <- c("(Intercept)", "groupB")

  dana_obj <- list(
    X = matrix(rnorm(40), nrow = 10, dimnames = list(paste0("s", 1:10), feat_ids)),
    sdata = data.frame(sample_id = paste0("s", 1:10), group = rep(c("A", "B"), each = 5)),
    formula_rhs = ~ group,
    fit = data.frame(
      feat_id = rep(feat_ids, each = 2),
      Coefficient = rep(coef_terms, times = 4),
      Estimate = rnorm(8),
      `Pr(>|t|)` = runif(8, max = 0.9), # Forces failure as no P-val > 0.95
      stringsAsFactors = FALSE
    ),
    lrt = data.frame(),
    ranef = data.frame()
  )
  class(dana_obj) <- "dana"

  expect_warning(
    adjust_pval(dana_obj, padj_method = "storey", verbose = FALSE),
    "qvalue failed; falling back to BH"
  )
})

test_that("adjust_pval adds LRT-adjusted p-values and maps them to dana_obj$fit", {
  skip_if_not_installed("qvalue")

  dana_obj <- list(
    X = matrix(rnorm(40), nrow = 10),
    sdata = data.frame(sample_id = paste0("s", 1:10), group = rep(c("A", "B"), each = 5)),
    formula_rhs = ~ group,
    fit = data.frame(
      feat_id = rep(c("f1", "f2"), each = 2),
      Coefficient = rep(c("(Intercept)", "groupB"), 2),
      Estimate = rnorm(4),
      `Pr(>|t|)` = runif(4),
      stringsAsFactors = FALSE
    ),
    lrt = data.frame(
      feat_id = c("f1", "f2"),
      term = "group",
      `Pr(>Chisq)` = runif(2)
    ),
    ranef = data.frame()
  )
  class(dana_obj) <- "dana"

  result <- adjust_pval(dana_obj, padj_method = "BH", padj_method_LRT = "BH", verbose = FALSE)

  expect_true("padj_BH_group_LRT" %in% colnames(result$fit))
  expect_true("padj_BH_group" %in% colnames(result$lrt))
  expect_true(all(is.numeric(result$fit$padj_BH_group_LRT), na.rm = TRUE))
})

test_that("adjust_pval supports interaction terms with colon (:) in LRT", {
  dana_obj <- list(
    X = matrix(rnorm(40), nrow = 10),
    sdata = data.frame(sample_id = paste0("s", 1:10), group = rep(c("A", "B"), each = 5)),
    formula_rhs = ~ group * time,
    fit = data.frame(
      feat_id = c("f1", "f1", "f2", "f2"),
      Coefficient = c("groupB", "groupB:timeT2", "groupB", "groupB:timeT2"),
      Estimate = rnorm(4),
      `Pr(>|t|)` = runif(4),
      stringsAsFactors = FALSE
    ),
    lrt = data.frame(
      feat_id = c("f1", "f2"),
      term = "group:time",
      `Pr(>Chisq)` = runif(2)
    ),
    ranef = data.frame()
  )
  class(dana_obj) <- "dana"

  result <- adjust_pval(dana_obj, padj_method = "BH", padj_method_LRT = "BH", verbose = FALSE)

  expect_true("padj_BH_group:time_LRT" %in% colnames(result$fit))
  expect_true(any(!is.na(result$fit$`padj_BH_group:time_LRT`)))
})

test_that("adjust_pval excludes random effect LRT terms from merging into fit", {
  skip_if_not_installed("stringr")

  # Simulate fit and LRT results
  fit_df <- data.frame(
    feat_id = rep(paste0("feat", 1:3), each = 2),
    Coefficient = rep(c("(Intercept)", "groupB"), 3),
    Estimate = rnorm(6),
    `Pr(>|t|)` = runif(6),
    stringsAsFactors = FALSE
  )

  lrt_df <- data.frame(
    feat_id = rep(paste0("feat", 1:3), 2),
    term = rep(c("group", "group|subject"), each = 3), # Mix fixed + random effect
    `Pr(>Chisq)` = runif(6),
    stringsAsFactors = FALSE
  )

  dana_obj <- list(
    fit = fit_df,
    lrt = lrt_df
  )
  class(dana_obj) <- "dana"

  # Run adjustment
  dana_out <- adjust_pval(dana_obj,
                          padj_by = "terms",
                          padj_method = "BH",
                          padj_method_LRT = "BH",
                          verbose = FALSE)

  # Check LRT adjustments exist for both terms in dana$lrt
  expect_true("padj_BH_group" %in% colnames(dana_out$lrt))
  expect_true("padj_BH_group|subject" %in% colnames(dana_out$lrt))

  # Check that only fixed-effect LRT term ("group") is merged into dana$fit
  expect_true("padj_BH_group_LRT" %in% colnames(dana_out$fit))
  expect_false("padj_BH_group|subject_LRT" %in% colnames(dana_out$fit))

  # Optional: check that padj values are numeric and NA only for unmatched
  expect_type(dana_out$fit$padj_BH_group_LRT, "double")
})

