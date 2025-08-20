test_that("build_phyloseq returns expected list structure with full input", {
  skip_if_not_installed("phyloseq")
  skip_if_not_installed("ape")

  mock_X <- matrix(
    c(10, 0, 5, 3, 1, 7),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("sample1", "sample2"), c("ASV1", "ASV2", "ASV3"))
  )

  mock_sample_data <- data.frame(
    sample_id = c("sample1", "sample2"),
    group = c("A", "B"),
    row.names = c("sample1", "sample2")
    )

  mock_taxa_table <- data.frame(
    Kingdom = c("Bacteria", "Bacteria", "Bacteria"),
    Genus = c("GenusA", "GenusB", "Unknown"),
    row.names = c("ASV1", "ASV2", "ASV3")
    )

  mock_tree <- ape::rtree(n = 3)
  mock_tree$tip.label <- rownames(mock_taxa_table)

  result <- build_phyloseq(
    X = mock_X,
    sample_data = mock_sample_data,
    taxa_table = mock_taxa_table,
    phylo_tree = mock_tree,
    taxa_in_rows = FALSE,
    verbose = FALSE
  )

  expect_type(result, "list")
  expect_true("asv" %in% names(result))
  expect_s4_class(result$asv, "phyloseq")
  expect_true("genus" %in% names(result))
  expect_s4_class(result$genus, "phyloseq")
})

test_that("build_phyloseq errors on mismatched data dimensions", {
  mock_X <- matrix(
    c(1, 2),
    nrow = 1,
    dimnames = list("sampleA", c("ASV1", "ASV2"))
  )

  mock_sample_data <- data.frame(
    sample_id = "sampleB",
    condition = "healthy",
    row.names = "sampleB",
    stringsAsFactors = FALSE
  )

  expect_error(
    build_phyloseq(
      X = mock_X,
      sample_data = mock_sample_data,
      taxa_table = data.frame(Kingdom = "Bacteria", row.names = "ASV1"),
      taxa_in_rows = FALSE,
      verbose = FALSE
    ),
    regexp = "Number of features"
  )
})

test_that("build_phyloseq errors on mismatched sample_id", {
  mock_X <- matrix(
    c(10, 0, 5, 3, 1, 7),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("sample1", "sample2"), c("ASV1", "ASV2", "ASV3"))
  )

  mock_sample_data <- data.frame(
    sample_id = c("sample3", "sample4"),
    group = c("A", "B"),
    row.names = c("sample3", "sample4"),
    stringsAsFactors = FALSE
  )

  mock_taxa_table <- data.frame(
    Kingdom = c("Bacteria", "Bacteria", "Bacteria"),
    Genus = c("GenusA", "GenusB", "Unknown"),
    row.names = c("ASV1", "ASV2", "ASV3"),
    stringsAsFactors = FALSE
  )

  mock_tree <- ape::rtree(n = 3)
  mock_tree$tip.label <- rownames(mock_taxa_table)

  expect_error(
    build_phyloseq(
      X = mock_X,
      sample_data = mock_sample_data,
      taxa_table = mock_taxa_table,
      phylo_tree = mock_tree,
      taxa_in_rows = FALSE,
      verbose = FALSE
    ),
    regexp = "set as row names"
  )
})


test_that("build_phyloseq warns about long feature labels", {
  skip_if_not_installed("phyloseq")

  long_ids <- c("ACGTGCTAGCAGTCA", "TGCAGATCGTAGGCA", "CGTAGCTAGCTAGCT")
  mock_X <- matrix(
    c(10, 20, 30, 40, 50, 60),
    nrow = 2,
    dimnames = list(c("sample1", "sample2"), long_ids)
  )

  mock_sample_data <- data.frame(
    sample_id = c("sample1", "sample2"),
    group = c("X", "Y"),
    row.names = c("sample1", "sample2")
    )

  mock_taxa_table <- data.frame(
    Kingdom = rep("Bacteria", 3),
    Genus = c("GenusA", "GenusB", "GenusC"),
    row.names = long_ids
    )

  mock_tree <- ape::rtree(n = 3)
  mock_tree$tip.label <- rownames(mock_taxa_table)

  expect_warning(
    result <- build_phyloseq(
      X = mock_X,
      sample_data = mock_sample_data,
      taxa_table = mock_taxa_table,
      phylo_tree = mock_tree,
      taxa_in_rows = FALSE,
      verbose = TRUE
    ),
    regexp = "names are too long"
  )

  expect_s4_class(result$asv, "phyloseq")
})

test_that("build_phyloseq filters out taxa with zero abundance", {
  skip_if_not_installed("phyloseq")

  mock_X <- matrix(
    c(0, 0, 10, 2, 5, 3),
    nrow = 2,
    dimnames = list(c("sample1", "sample2"), c("ASV1", "ASV2", "ASV3"))
  )

  mock_sample_data <- data.frame(
    sample_id = c("sample1", "sample2"),
    group = c("A", "B"),
    row.names = c("sample1", "sample2"),
    stringsAsFactors = FALSE
  )

  mock_taxa_table <- data.frame(
    Kingdom = rep("Bacteria", 3),
    Genus = c("GenusA", "GenusB", "GenusC"),
    row.names = c("ASV1", "ASV2", "ASV3"),
    stringsAsFactors = FALSE
  )

  mock_tree <- ape::rtree(n = 3)
  mock_tree$tip.label <- rownames(mock_taxa_table)


  result <- build_phyloseq(
    X = mock_X,
    sample_data = mock_sample_data,
    taxa_table = mock_taxa_table,
    phylo_tree = mock_tree,
    taxa_in_rows = FALSE,
    verbose = FALSE
  )

  retained_taxa <- phyloseq::taxa_names(result$asv)
  expect_false("ASV1" %in% retained_taxa)
  expect_true("ASV2" %in% retained_taxa)
})

test_that("build_phyloseq filters out samples with zero abundance", {
  skip_if_not_installed("phyloseq")

  mock_X <- matrix(
    c(0, 0, 0, 5, 0, 3),
    nrow = 2,
    dimnames = list(c("sample1", "sample2"), c("ASV1", "ASV2", "ASV3"))
  )

  mock_sample_data <- data.frame(
    sample_id = c("sample1", "sample2"),
    group = c("A", "B"),
    row.names = c("sample1", "sample2"),
    stringsAsFactors = FALSE
  )

  mock_taxa_table <- data.frame(
    Kingdom = rep("Bacteria", 3),
    Genus = c("GenusA", "GenusB", "GenusC"),
    row.names = c("ASV1", "ASV2", "ASV3"),
    stringsAsFactors = FALSE
  )

  mock_tree <- ape::rtree(n = 3)
  mock_tree$tip.label <- rownames(mock_taxa_table)

  result <- build_phyloseq(
    X = mock_X,
    sample_data = mock_sample_data,
    taxa_table = mock_taxa_table,
    phylo_tree = mock_tree,
    taxa_in_rows = FALSE,
    verbose = FALSE
    )

  retained_samples <- phyloseq::sample_names(result$asv)
  expect_false("sample1" %in% retained_samples)
  expect_true("sample2" %in% retained_samples)
})
