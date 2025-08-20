test_that("process_ms runs with default parameters", {
  set.seed(1)
  X <- matrix(sample(c(0:10), size = 80, replace = TRUE),
              nrow = 20, ncol = 4,
              dimnames = list(paste0("sample", 1:20), paste0("feat", 1:4)))

  expect_warning(
    result <- process_ms(X, verbose = FALSE),
    regexp = "contains NAs"
  )

  expect_type(result, "list")
  expect_named(result, c("X_names", "X_processed"))
  expect_s3_class(result$X_names, "data.frame")
  expect_true(all(grepl("^feat_", colnames(result$X_processed))))
})


test_that("process_ms removes samples by pattern", {
  X <- matrix(abs(rnorm(80)), nrow = 20,
              dimnames = list(c(paste0("qc_", 1:5), paste0("sample", 1:15)),
                              paste0("feat", 1:4)))

  expect_warning(
    result <- process_ms(X, remove_ids = "^qc_", verbose = FALSE),
    regexp = "contains NAs"
  )

  expect_equal(nrow(result$X_processed), 15)
  expect_false(any(grepl("^qc_", rownames(result$X_processed))))
})


test_that("process_ms performs log transformation", {
  X <- matrix(2^sample(1:10, 80, replace = TRUE),
              nrow = 20, dimnames = list(NULL, paste0("feat", 1:4)))

  rownames(X) <- paste0("sample", seq_len(nrow(X)))

  expect_warning(
    result <- process_ms(X, transform = "log", log_base_num = 2, verbose = FALSE),
    regexp = "contains NAs"
  )

  expect_equal(result$X_processed[1, 1], log2(X[1, 1]))
})


test_that("process_ms applies min_val imputation", {
  set.seed(165)
  X <- matrix(sample(c(1:10, NA), size = 80, replace = TRUE,
                     prob = c(rep(0.08, 10), 0.2)),
              nrow = 20)

  colnames(X) <- paste0("feat_", seq_len(ncol(X)))
  rownames(X) <- paste0("sample", seq_len(nrow(X)))

  result <- process_ms(X, impute = "min_val", min_val_factor = 2, verbose = TRUE)

  expect_false(any(is.na(result$X_processed)))
  expect_equal(dim(result$X_processed), c(20, 2))
})

test_that("process_ms performs QRILC imputation with transform override", {
  skip_if_not_installed("imputeLCMD")

  set.seed(165)
  X <- matrix(sample(c(5:15, NA), size = 80, replace = TRUE,
                     prob = c(rep(0.08, 11), 0.12)),
              nrow = 20)

  colnames(X) <- paste0("feat_", seq_len(ncol(X)))
  rownames(X) <- paste0("sample", seq_len(nrow(X)))

  expect_warning(
    result <- process_ms(X, transform = "none", impute = "QRILC", verbose = TRUE),
    regexp = "Overriding.*transform.*log"
  )

  expect_false(any(is.na(result$X_processed)))
  expect_equal(dim(result$X_processed), c(20, 4))
})


test_that("process_ms validates feature renaming", {
  X <- matrix(runif(80), nrow = 20,
              dimnames = list(NULL, paste0("M", 1:4)))

  colnames(X) <- paste0("complicatedid_", seq_len(ncol(X)))
  rownames(X) <- paste0("sample", seq_len(nrow(X)))

  expect_warning(
    result <- process_ms(X, impute = "none", verbose = FALSE),
    regexp = "contains NAs"
  )

  expect_true(all(grepl("^feat_", colnames(result$X_processed))))
  expect_equal(ncol(result$X_processed), nrow(result$X_names))
})


test_that("process_ms validates input arguments", {
  X <- matrix(runif(80), nrow = 20)

  colnames(X) <- paste0("feat_", seq_len(ncol(X)))
  rownames(X) <- paste0("sample", seq_len(nrow(X)))

  expect_error(process_ms(X, transform = "log", log_base_num = NA),
               regexp = "log_base_num.*must be a numeric")

  expect_error(process_ms(X, impute = "min_val", min_val_factor = 0.5),
               regexp = "min_val_factor.*>= 1")

  expect_error(process_ms(X, seed = "xyz"),
               regexp = "seed.*numeric value")
})


test_that("process_ms warns if imputation skipped due to no NAs", {
  X <- matrix(1:80, nrow = 20)

  colnames(X) <- paste0("feat_", seq_len(ncol(X)))
  rownames(X) <- paste0("sample", seq_len(nrow(X)))

  expect_message(
    result <- process_ms(X, impute = "min_val", verbose = TRUE),
    regexp = "no missing values"
  )

  expect_equal(dim(result$X_processed), c(20, 4))
})
