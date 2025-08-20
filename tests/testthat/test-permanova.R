test_that("permanova runs joint and independent PERMANOVA correctly", {
  skip_if_not_installed("vegan")
  skip_if_not_installed("permute")

  # Mock data
  X <- matrix(rnorm(40), nrow = 10,
              dimnames = list(paste0("sample", 1:10), paste0("feat", 1:4)))
  sample_data <- data.frame(
    sample_id = rownames(X),
    group = rep(c("A", "B"), each = 5),
    age = rep(20:29, length.out = 10),
    row.names = rownames(X),
    stringsAsFactors = FALSE
  )

  # Simple control structures
  dist_control <- list(method = "euclidean")
  perm_control <- list(
    joint_terms = list(control = permute::how(blocks = NULL, nperm = 9)),
    group = list(control = permute::how(blocks = NULL, nperm = 9)),
    age = list(control = permute::how(blocks = NULL, nperm = 9))
  )

  result <- permanova(
    X = X,
    sample_data = sample_data,
    formula_rhs = ~ group + age,
    dist_control = dist_control,
    perm_control = perm_control,
    independent = TRUE,
    platform = "ms",
    assay = "lipidomics",
    seed = 42,
    verbose = FALSE
  )

  expect_type(result, "list")
  expect_named(result, c("X_dist", "perm_matrix_joint", "permanova_joint", "permanova_indep"))
  expect_s3_class(result$X_dist, "dist")
  expect_s3_class(result$permanova_joint, "anova")
  expect_s3_class(result$permanova_indep, "data.frame")
  expect_true(all(c("platform", "assay") %in% colnames(result$permanova_indep)))
})


test_that("permanova returns consistent results with seed", {
  skip_if_not_installed("vegan")
  skip_if_not_installed("permute")

  X <- matrix(rnorm(40), nrow = 10)
  rownames(X) <- paste0("sample", 1:10)

  sample_data <- data.frame(
    sample_id = rownames(X),
    group = rep(c("A", "B"), each = 5),
    row.names = rownames(X)
  )

  perm_control <- list(
    joint_terms = list(control = permute::how(blocks = NULL, nperm = 9))
  )

  expect_warning(
    r1 <- permanova(X, sample_data, ~ group,
                    perm_control = perm_control,
                    seed = 123, platform = "nmr"),
    regexp = "independent = TRUE"
  )

  expect_warning(
    r2 <- permanova(X, sample_data, ~ group,
                    perm_control = perm_control,
                    seed = 123, platform = "nmr"),
    regexp = "independent = TRUE"
  )

  expect_equal(r1$permanova_joint$`Pr(>F)`, r2$permanova_joint$`Pr(>F)`)
})


test_that("permanova errors with formula that includes response", {
  skip_if_not_installed("vegan")

  X <- matrix(rnorm(40), nrow = 10)
  rownames(X) <- paste0("sample", 1:10)
  sample_data <- data.frame(group = rep(c("A", "B"), each = 5),
                            row.names = rownames(X))

  expect_error(
    permanova(X, sample_data, group ~ 1, platform = "ms"),
    regexp = "one-sided formula"
  )
})


test_that("permanova warns when specific perm_control is missing but independent = TRUE", {
  skip_if_not_installed("vegan")
  skip_if_not_installed("permute")

  X <- matrix(rnorm(40), nrow = 10)
  rownames(X) <- paste0("sample", 1:10)
  sample_data <- data.frame(
    sample_id = rownames(X),
    group = rep(c("A", "B"), each = 5),
    row.names = rownames(X)
  )

  perm_control <- list(
    joint_terms = list(control = permute::how(blocks = NULL, nperm = 9))
  )

  expect_warning(
    permanova(X, sample_data, ~ group, perm_control = perm_control,
              independent = TRUE, platform = "ngs"),
    regexp = "independent = TRUE"
  )
})


test_that("permanova supports vegan::vegdist methods", {
  skip_if_not_installed("vegan")

  X <- matrix(runif(50), nrow = 10)
  rownames(X) <- paste0("sample", 1:10)

  sample_data <- data.frame(
    sample_id = rownames(X),
    group = factor(rep(c("A", "B"), each = 5)),
    row.names = rownames(X)
  )

  dist_control <- list(method = "bray")  # vegan::vegdist
  perm_control <- list(joint_terms = list(control = permute::how(blocks = NULL, nperm = 9)))

  expect_warning(
    result <- permanova(X, sample_data, ~ group,
                        dist_control = dist_control,
                        perm_control = perm_control,
                        platform = "ms",
                        verbose = FALSE),
    regexp = "independent = TRUE"
  )

  expect_s3_class(result$X_dist, "dist")
})

test_that("permanova handles NA in predictor variables when na.action = na.omit", {
  set.seed(123)

  # Fake data
  X <- matrix(rnorm(50 * 5), nrow = 50)
  rownames(X) <- paste0("sample", seq_len(50))

  sample_data <- data.frame(
    sample_id = rownames(X),
    group = sample(c("A", "B"), 50, replace = TRUE),
    age = rnorm(50),
    condition = sample(c("healthy", "disease"), 50, replace = TRUE),
    stringsAsFactors = TRUE
  )
  rownames(sample_data) <- sample_data$sample_id

  # Inject NAs into one variable
  sample_data$age[c(3, 12, 45)] <- NA

  # Setup inputs
  formula_rhs <- ~ group + age + condition
  dist_control <- list(method = "euclidean")
  perm_control <- list(
    joint_terms = list(control = permute::how(blocks = NULL, nperm = 99))
  )

  expect_warning(
    expect_warning(
      result <- permanova(X, sample_data, formula_rhs,
                dist_control = dist_control,
                perm_control = perm_control,
                independent = TRUE,
                platform = "nmr",
                assay = "test",
                seed = 42, na.action = na.omit,
                verbose = FALSE),
      regexp = "perm_control"),
    regexp = "Missing data")

  expect_type(result, "list")
  expect_named(result, c("X_dist", "perm_matrix_joint", "permanova_joint", "permanova_indep"))
  expect_s3_class(result$X_dist, "dist")
  expect_s3_class(result$permanova_joint, "anova")
  expect_s3_class(result$permanova_indep, "data.frame")

  # Check that samples with NA in age were dropped
  expect_lt(attr(result$X_dist, "Size"), 50)  # X_dist reflects number of samples after NA removal
})
