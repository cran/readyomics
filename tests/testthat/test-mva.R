test_that("mva runs PCA and returns expected output structure", {
  skip_if_not_installed("ropls")
  skip_if_not_installed("ggplot2")

  set.seed(123)
  mock_X <- matrix(
    rnorm(40),
    nrow = 10,
    dimnames = list(paste0("sample", 1:10), paste0("feat", 1:4))
  )

  sample_data <- data.frame(
    sample_id = rownames(mock_X),
    group = rep(c("A", "B"), each = 5),
    batch = rep(1:2, times = 5),
    row.names = rownames(mock_X),
    stringsAsFactors = FALSE
  )

  result <- mva(
    X = mock_X,
    sample_data = sample_data,
    group_colour = "group",
    group_shape = "batch",
    plot_title = "Test PCA Plot",
    predI = 2,  # PCA: set components
    verbose = FALSE
  )

  expect_type(result, "list")
  expect_named(result, c("ropls_obj", "scores_plot"))
  expect_s4_class(result$ropls_obj, "opls")
  expect_s3_class(result$scores_plot, "ggplot")
})


test_that("mva errors when rownames do not match between X and sample_data", {
  skip_if_not_installed("ropls")

  X <- matrix(rnorm(30), nrow = 5)
  rownames(X) <- paste0("s", 1:5)

  sample_data <- data.frame(
    sample_id = paste0("x", 1:5),  # mismatch
    group = rep(c("A", "B"), length.out = 5),
    row.names = paste0("x", 1:5)
  )

  expect_error(
    mva(X = X, sample_data = sample_data, predI = 2),
    regexp = "sample_id.*row names|match"
  )
})


test_that("mva handles NULL group_colour and group_shape", {
  skip_if_not_installed("ropls")
  skip_if_not_installed("ggplot2")

  X <- matrix(rnorm(40), nrow = 10)
  rownames(X) <- paste0("sample", 1:10)

  sample_data <- data.frame(
    sample_id = rownames(X),
    row.names = rownames(X)
  )

  res <- mva(
    X = X,
    sample_data = sample_data,
    predI = 2,
    verbose = FALSE
  )

  expect_s3_class(res$scores_plot, "ggplot")
})


test_that("mva throws error when invalid colour or shape column is given", {
  skip_if_not_installed("ropls")

  X <- matrix(rnorm(40), nrow = 10)
  rownames(X) <- paste0("sample", 1:10)

  sample_data <- data.frame(
    sample_id = rownames(X),
    group = rep(c("A", "B"), each = 5),
    row.names = rownames(X)
  )

  expect_error(
    mva(
      X = X,
      sample_data = sample_data,
      group_colour = "nonexistent",
      predI = 2
    ),
    regexp = "Cannot find"
  )
})


test_that("mva returns correct axis labels and title in ggplot object", {
  skip_if_not_installed("ropls")
  skip_if_not_installed("ggplot2")

  X <- matrix(rnorm(40), nrow = 10)
  rownames(X) <- paste0("sample", 1:10)

  sample_data <- data.frame(
    sample_id = rownames(X),
    group = rep(c("A", "B"), each = 5),
    row.names = rownames(X)
  )

  plot_title <- "My PCA Plot"

  result <- mva(
    X = X,
    sample_data = sample_data,
    group_colour = "group",
    plot_title = plot_title,
    predI = 2
  )

  p <- result$scores_plot
  expect_s3_class(p, "ggplot")
  expect_match(p$labels$title, plot_title)
  expect_true(is.call(p$labels$x))
  expect_true(is.call(p$labels$y))
})
