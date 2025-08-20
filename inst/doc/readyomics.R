## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval = FALSE-------------------------------------------------------------
# install.packages("readyomics")

## ----setup--------------------------------------------------------------------
library(readyomics)

## -----------------------------------------------------------------------------
# Raw matrix of ASV counts
asv_counts <- read.csv(system.file("extdata", "asv_raw_counts.csv", package = "readyomics"), 
                       check.names = FALSE, 
                       row.names = 1)

# Taxonomy table
taxa <- read.csv(system.file("extdata", "taxonomy.csv", package = "readyomics"), 
                 check.names = FALSE, 
                 row.names = 1)

# Sample data
sample_data <- read.csv(system.file("extdata", "sample_data.csv", package = "readyomics"), 
                        check.names = FALSE, 
                        row.names = 1)

## ----eval = FALSE-------------------------------------------------------------
# head(asv_counts)

## ----echo = FALSE-------------------------------------------------------------
knitr::kable(head(asv_counts))

## ----eval = FALSE-------------------------------------------------------------
# head(taxa)

## ----echo = FALSE-------------------------------------------------------------
knitr::kable(head(taxa))

## ----eval = FALSE-------------------------------------------------------------
# head(sample_data)

## ----echo = FALSE-------------------------------------------------------------
knitr::kable(head(sample_data))

## -----------------------------------------------------------------------------
sample_data <- sample_data |>
  dplyr::mutate(sample_id = rownames(sample_data),
         groups_liver = factor(groups_liver, levels = c("FL", "FL_HS")),
         sex = factor(sex),
         PPI = factor(PPI),
         smoking = factor(smoking),
         alcohol_imp = factor(alcohol_imp))

## -----------------------------------------------------------------------------
# Samples as rows required for process_ngs
asv_counts <- t(asv_counts)

# Process asv_counts data
asv_ready <- process_ngs(X = asv_counts, 
                         sample_data = sample_data, 
                         taxa_table = taxa, 
                         normalise = "none",
                         transform = "clr", 
                         eco_phyloseq = FALSE)

## -----------------------------------------------------------------------------
pca <- mva(X = asv_ready$X_processed, 
           sample_data = asv_ready$sdata_final, 
           group_colour = "groups_liver", 
           plot_title = "Beta diversity (Aitchison)")

## -----------------------------------------------------------------------------
pca$scores_plot

## -----------------------------------------------------------------------------
# Define model
rhs_model <- ~ groups_liver + alcohol_imp + sex + age + PPI + smoking

# Run PERMANOVA
pova <- permanova(X = asv_ready$X_processed, 
                  sample_data = asv_ready$sdata_final,
                  formula_rhs = rhs_model, 
                  platform = "ngs", 
                  assay = "Taxonomy", 
                  seed = 165)

## ----eval = FALSE-------------------------------------------------------------
# pova$permanova_joint

## ----echo = FALSE-------------------------------------------------------------
knitr::kable(pova$permanova_joint)

## ----eval = FALSE-------------------------------------------------------------
# pova$permanova_indep

## ----echo = FALSE-------------------------------------------------------------
knitr::kable(pova$permanova_indep)

## -----------------------------------------------------------------------------
# future::plan(multisession, workers = 4) # For running 4 processes in parallel
# progressr::handlers(global = TRUE) # To show progress bar
dana_asv <- dana(X = asv_ready$X_processed, 
                 sample_data = asv_ready$sdata_final,
                 formula_rhs = rhs_model, 
                 term_LRT = "groups_liver",
                 platform = "ngs")
# plan(sequential) # To disable parallel processing

## ----eval = FALSE-------------------------------------------------------------
# head(dana_asv$fit)

## ----echo = FALSE-------------------------------------------------------------
knitr::kable(head(dana_asv$fit))

## ----eval = FALSE-------------------------------------------------------------
# head(dana_asv$lrt)

## ----echo = FALSE-------------------------------------------------------------
knitr::kable(head(dana_asv$lrt))

## -----------------------------------------------------------------------------
dana_asv <- dana_asv |>
  adjust_pval(padj_by = "terms", 
              padj_method = "BH", 
              padj_method_LRT = "BH")

## ----eval = FALSE-------------------------------------------------------------
# head(dana_asv$fit)

## ----echo = FALSE-------------------------------------------------------------
knitr::kable(head(dana_asv$fit))

## ----eval = FALSE-------------------------------------------------------------
# head(dana_asv$lrt)

## ----echo = FALSE-------------------------------------------------------------
knitr::kable(head(dana_asv$lrt))

## -----------------------------------------------------------------------------
dana_asv <- dana_asv |>
  add_taxa(taxa_table = taxa, 
           taxa_rank = "asv")

## ----eval = FALSE-------------------------------------------------------------
# head(dana_asv$fit)

## ----echo = FALSE-------------------------------------------------------------
knitr::kable(head(dana_asv$fit))

## ----error = TRUE-------------------------------------------------------------
try({
# Generates error due to lack of significant features:
dana_asv <- dana_asv |>
  ready_plots(term_name = "groups_liver", # Formula fit term of interest
              pval_match = "groups_liver_LRT", # LRT adjusted P values will be used
              sdata_var = "groups_liver", # Grouping variable for individual feature plots
              group_colours = c(FL = "#4daf4a", # Optional color customization
                                FL_HS = "#377eb8"))
})

## -----------------------------------------------------------------------------
dana_asv <- dana_asv |>
  ready_plots(term_name = "groups_liver",
              pval_match = "Pr", # Nominal P value for illustration purposes
              sdata_var = "groups_liver", 
              group_colours = c(FL = "#4daf4a", 
                                FL_HS = "#377eb8"))

## -----------------------------------------------------------------------------
dana_asv$plots

## -----------------------------------------------------------------------------
# Compute alpha diversity measures for the ASV phyloseq object
alpha <- phyloseq::estimate_richness(asv_ready$phyloseq_raw$asv) |>
  dplyr::select(-matches("ace|chao|fisher")) |>
  scale()

# Calculate alpha diversity differences among groups with dana()
dana_alpha <- dana(X = alpha, 
                   sample_data = asv_ready$sdata_final,
                   formula_rhs = rhs_model, 
                   term_LRT = "groups_liver",
                   platform = "ngs")

## -----------------------------------------------------------------------------
# CLR-transform genus rank counts table
genus_ready <- process_ngs(X = as.data.frame(asv_ready$phyloseq_raw$genus@otu_table), 
                           sample_data = asv_ready$sdata_final, 
                           taxa_table = taxa, 
                           normalise = "none",
                           transform = "clr",
                           raw_phyloseq = FALSE,
                           eco_phyloseq = FALSE,
                           verbose = FALSE)

# Run dana
dana_genus <- dana(X = genus_ready$X_processed, 
                   sample_data = genus_ready$sdata_final,
                   formula_rhs = rhs_model, 
                   term_LRT = "groups_liver",
                   platform = "ngs") |>
  adjust_pval(padj_by = "terms", 
              padj_method = "BH", 
              padj_method_LRT = "BH") |>
  add_taxa(taxa_table = taxa, 
           taxa_rank = "genus") |>
  ready_plots(term_name = "groups_liver",
              pval_match = "Pr", # Nominal P value for illustration purposes
              sdata_var = "groups_liver",
              group_colours = c(FL = "#4daf4a",
                                FL_HS = "#377eb8"))

# Inspect plots
dana_genus$plots

