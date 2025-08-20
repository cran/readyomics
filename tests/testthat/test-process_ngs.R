test_that("process_ngs runs with minimal valid input", {
  skip_if_not_installed("phyloseq")
  skip_if_not_installed("zCompositions")

  # Create 5x5 count matrix
  set.seed(123)
  mock_X <- matrix(
    sample(0:1000, 25, replace = TRUE),
    nrow = 5,
    dimnames = list(paste0("sample", 1:5), paste0("ASV", 1:5))
  )

  mock_sample_data <- data.frame(
    sample_id = paste0("sample", 1:5),
    load = c(1e5, 2e5, 1e4, 5e4, 1.5e5),
    condition = rep(c("A", "B"), length.out = 5),
    row.names = paste0("sample", 1:5)
  )

  mock_taxa_table <- data.frame(
    Kingdom = rep("Bacteria", 5),
    Genus = paste0("Genus", 1:5),
    row.names = paste0("ASV", 1:5)
  )

  mock_tree <- ape::rtree(5)
  mock_tree$tip.label <- rownames(mock_taxa_table)

  result <- process_ngs(
    X = mock_X,
    sample_data = mock_sample_data,
    taxa_table = mock_taxa_table,
    phylo_tree = mock_tree,
    normalise = "load",
    load_colname = "load",
    transform = "none",
    verbose = FALSE
    )

  expect_type(result, "list")
  expect_true("X_processed" %in% names(result))
  expect_true("phyloseq_raw" %in% names(result))
  expect_s4_class(result$phyloseq_raw$asv, "phyloseq")
  expect_s4_class(result$phyloseq_eco$asv, "phyloseq")
  expect_true(nrow(result$sdata_final) > 0)
})

test_that("process_ngs throws error with invalid load column", {
  mock_X <- matrix(
    c(10, 20, 30, 40, 50, 60),
    nrow = 2,
    dimnames = list(c("sample1", "sample2"), c("ASV1", "ASV2", "ASV3"))
  )
  mock_sample_data <- data.frame(
    sample_id = c("sample1", "sample2"),
    group = c("X", "Y"),
    row.names = c("sample1", "sample2")
  )

  expect_error(
    process_ngs(
      X = mock_X,
      sample_data = mock_sample_data,
      normalise = "load",
      load_colname = "missing_column"
    ),
    "not found in 'sample_data'"
  )
})

test_that("process_ngs warns if microbial load is below threshold", {
  skip_if_not_installed("zCompositions")

  set.seed(123)

  n_samples <- 20
  n_features <- 50

  # Build 20x50 matrix: each column with 0–20 zeros
  mock_X <- matrix(NA_integer_, nrow = n_samples, ncol = n_features)
  for (j in seq_len(n_features)) {
    n_zeros <- sample(0:20, 1)
    non_zeros <- sample(1:1000, n_samples - n_zeros, replace = TRUE)
    col_vals <- c(rep(0L, n_zeros), non_zeros)
    mock_X[, j] <- sample(col_vals)
  }

  rownames(mock_X) <- paste0("sample", seq_len(n_samples))
  colnames(mock_X) <- paste0("ASV", seq_len(n_features))

  mock_sample_data <- data.frame(
    sample_id = paste0("sample", seq_len(n_samples)),
    load = sample(100:100000, size = n_samples),
    row.names = paste0("sample", seq_len(n_samples))
  )

  expect_warning(
    process_ngs(
      X = mock_X,
      sample_data = mock_sample_data,
      normalise = "load",
      load_colname = "load",
      transform = "log",
      raw_phyloseq = FALSE,
      eco_phyloseq = FALSE,
      verbose = FALSE
    ),
    "Microbial load for some samples is <"
  )
})

test_that("process_ngs returns intermediates when return_all = TRUE", {
  skip_if_not_installed("zCompositions")

  set.seed(123)

  n_samples <- 20
  n_features <- 50

  # Build 20x50 matrix: each column with 0–20 zeros
  mock_X <- matrix(NA_integer_, nrow = n_samples, ncol = n_features)
  for (j in seq_len(n_features)) {
    n_zeros <- sample(0:20, 1)
    non_zeros <- sample(1:1000, n_samples - n_zeros, replace = TRUE)
    col_vals <- c(rep(0L, n_zeros), non_zeros)
    mock_X[, j] <- sample(col_vals)
  }

  rownames(mock_X) <- paste0("sample", seq_len(n_samples))
  colnames(mock_X) <- paste0("ASV", seq_len(n_features))

  mock_sample_data <- data.frame(
    sample_id = paste0("sample", seq_len(n_samples)),
    load = sample(10000:100000, size = n_samples),
    row.names = paste0("sample", seq_len(n_samples))
  )

  result <- process_ngs(
      X = mock_X,
      sample_data = mock_sample_data,
      normalise = "load",
      load_colname = "load",
      transform = "clr",
      raw_phyloseq = FALSE,
      eco_phyloseq = FALSE,
      return_all = TRUE,
      verbose = FALSE
    )


  expect_true("X_matched" %in% names(result))
  expect_true("X_norm" %in% names(result))
  expect_true("X_prev" %in% names(result))
})

test_that("process_ngs errors if all samples filtered by min_reads", {
  mock_X <- matrix(
    0,
    nrow = 2, ncol = 5,
    dimnames = list(c("sample1", "sample2"), paste0("ASV", 1:5))
  )

  mock_sample_data <- data.frame(
    sample_id = c("sample1", "sample2"),
    load = c(1e5, 1e5),
    row.names = c("sample1", "sample2")
  )

  expect_error(
    process_ngs(
      X = mock_X,
      sample_data = mock_sample_data,
      normalise = "load",
      load_colname = "load",
      raw_phyloseq = FALSE,
      eco_phyloseq = FALSE,
      verbose = FALSE
    ),
    "No samples remain after row names and min_reads filtering"
  )
})

