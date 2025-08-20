test_that("dana runs with basic linear model", {
  skip_if_not_installed("lmerTest")
  skip_if_not_installed("future")
  skip_if_not_installed("progressr")

  set.seed(123)
  mock_X <- matrix(
    rnorm(20), nrow = 5,
    dimnames = list(paste0("sample", 1:5), paste0("feat", 1:4))
  )

  mock_sample_data <- data.frame(
    sample_id = rownames(mock_X),
    group = factor(c("A", "A", "B", "B", "A")),
    row.names = rownames(mock_X)
  )

  expect_warning(
    result <- dana(
      X = mock_X,
      sample_data = mock_sample_data,
      formula_rhs = ~ group,
      platform = "ms",
      verbose = FALSE
      ),
    regexp = "Likelihood ratio test"
  )

  expect_s3_class(result, "dana")
  expect_true(all(c("X", "sdata", "fit", "ranef", "lrt", "errors") %in% names(result)))
  expect_true(nrow(result$fit) > 0)
  expect_true(nrow(result$ranef) == 0) # no random effects
  expect_true(nrow(result$lrt) == 0) # no LRT term
  expect_true(nrow(result$errors) == 0) # no model failures
})


test_that("dana handles LRT terms correctly", {
  skip_if_not_installed("lmerTest")
  skip_if_not_installed("future")

  mock_X <- matrix(rnorm(20), nrow = 5)
  rownames(mock_X) <- paste0("id", 1:5)
  colnames(mock_X) <- paste0("feat", 1:4)

  sample_data <- data.frame(
    sample_id = rownames(mock_X),
    group = factor(c("A", "A", "B", "B", "B")),
    subject = factor(c(1, 1, 2, 2, 3)),
    row.names = rownames(mock_X)
  )

  res <- dana(
    X = mock_X,
    sample_data = sample_data,
    formula_rhs = ~ group + (1 | subject),
    term_LRT = "group",
    platform = "nmr",
    verbose = FALSE
    )

  expect_s3_class(res, "dana")
  expect_true(nrow(res$fit) > 0)
  expect_true(nrow(res$lrt) > 0)
  expect_true("term" %in% colnames(res$lrt))
})


test_that("dana errors on bad formula_rhs", {
  mock_X <- matrix(rnorm(12), nrow = 4)

  sample_data <- data.frame(
    sample_id = paste0("s", 1:4),
    group = c("A", "A", "B", "B"),
    row.names = paste0("s", 1:4)
  )

  expect_error(
      dana(X = mock_X,
           sample_data = sample_data,
           formula_rhs = y ~ group,
           platform = "ms"),
      regexp = "one-sided formula"
    )
})


test_that("dana errors on unmatched formula terms", {
  mock_X <- matrix(rnorm(12), nrow = 4)
  sample_data <- data.frame(
    sample_id = paste0("s", 1:4),
    group = c("A", "A", "B", "B"),
    row.names = paste0("s", 1:4)
  )

  expect_error(
    dana(X = mock_X, sample_data = sample_data, formula_rhs = ~ dose, platform = "ngs"),
    regexp = "formula_rhs.*sample_data"
    )
})


test_that("dana errors if term_LRT not in formula_rhs", {
  mock_X <- matrix(rnorm(12), nrow = 4)
  sample_data <- data.frame(
    sample_id = paste0("s", 1:4),
    group = c("A", "A", "B", "B"),
    row.names = paste0("s", 1:4)
  )

  expect_error(
    dana(
        X = mock_X,
        sample_data = sample_data,
        formula_rhs = ~ group,
        term_LRT = "dose",
        platform = "ms"
      ),
      regexp = "term_LRT.*formula_rhs"
  )
})


test_that("dana handles model errors per feature", {
  mock_X <- matrix(
    c(
      10,  5,  8,  Inf,
      12,  6,  7,  4,
      11,  7,  9,  8,
      13,  5,  6,  1,
      12,  4,  7,  10
    ),
    nrow = 5,
    dimnames = list(paste0("sample", 1:5), paste0("feat", 1:4))
  )

  sample_data <- data.frame(
    sample_id = rownames(mock_X),
    group = factor(c("A", "A", "B", "B", "B")),
    row.names = rownames(mock_X)
  )

  expect_warning(
    expect_warning(
      result <- dana(
          X = mock_X,
          sample_data = sample_data,
          formula_rhs = ~ group,
          platform = "ms",
          verbose = FALSE
        ),
        regexp = "Likelihood ratio test"
    ),
    regexp = "Model failed for feature"
  )

  expect_s3_class(result, "dana")
  expect_true(nrow(result$errors) > 0)
  expect_true("error" %in% colnames(result$errors))
  expect_equal(result$errors$feat_id, "feat1")
})
