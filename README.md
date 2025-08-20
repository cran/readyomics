
# readyomics <img src="man/figures/logo.png" align="right" height="130" alt="" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/lmartinezgili/readyomics/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/lmartinezgili/readyomics/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

readyomics provides a pipeline for formatting, analyzing, and visualizing omics 
data - regardless of omics type (e.g. transcriptomics, proteomics, metabolomics, 
metagenomics).

It is designed for flexibility, reproducibility, and scalability across a wide 
range of study designs, with modular components for statistical analysis 
and visualization.

It includes tools to:

- Process data into analysis-ready.
- Perform multivariate analysis.
- Fit linear or mixed-effects models.
- Produce publication-quality plots.

## Essential standard terms used in the package
- **Platform**: the technology or instrument used to generate omics data, such as 
  next-generation sequencing (NGS), mass spectrometry (MS), or nuclear magnetic 
  resonance spectroscopy (NMR).
- **Feature(s)**: a general term for a biological variable that has been profiled
  with an omics platform, such as metabolites, lipids, genes, proteins, or 
  microbial taxa, depending on the assay.
- **Sample data** (!= *metadata*\*): biological or demographic information collected 
  for each study sample (e.g., experimental group, age, sex, BMI).

\**Note*: in its strict sense, metadata ("data about data") refers to information 
describing the context, structure, or properties of a dataset — such as acquisition date, 
instrument settings, plate/well ID, or run order. It does not refer to biological 
or demographic variables.
To avoid ambiguity, readyomics adopts the same convention as phyloseq, using the 
term *sample data* for variables describing the study samples.

## Main functionalities

### Data processing (normalisation, transformation, filtering)
- `process_ngs()`: process next-generation sequencing data.
- `process_ms()`: process MS or NMR data.
- `build_phyloseq()`: build phyloseq objects for metataxonomic data.

### Multivariate analysis
- `mva()`: PCA, PLS and OPLS models.
- `permanova()`: wrapper for `vegan::adonis2()` function with additional options 
  and summary results.

### Differential [abundance/expression] analysis
- `dana()`: fit feature-wise linear fixed or mixed effects models.
- `adjust_pval()`: methods to adjust nominal P-values on `dana()` result.
- `ready_plots()`: visualize `dana` results and significant features.

## Installation
``` r
install.packages("readyomics")
```

You can install the development version of readyomics from GitHub:

``` r
devtools::install_github("lmartinezgili/readyomics")
```

## Get started
### Types of input data
readyomics is as omics-agnostic and inclusive as possible. 

Raw omics files (e.g., .fastq, .mzML) must first be pre-processed through 
external commercial or open-source pipelines into a data matrix where each row 
corresponds to a sample and each column corresponds to a measured omics feature.

### Required files and format
- `X`: a .csv or .RDS table of omics data (samples in rows and features in columns).
- `sample_data` (or `sdata`): a .csv or .RDS table of study sample information (samples in rows).
- `sample_id` must be a column in `sample_data` and have unique ids for each sample.
- Row names in `X` and `sample_data` must match `sample_id` values, though order 
  can differ — readyomics functions will check and align automatically.

### Documentation and Examples
For tutorials, examples, and reference documentation, visit readyomics [website](https://lmartinezgili.github.io/readyomics/).
