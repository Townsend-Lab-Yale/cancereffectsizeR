# Introduction
CancereffectsizeR is an R package that estimates the selection intensity of somatic mutations in cancer from collections of tumor sequencing data. In brief, each tumor sample's trinucleotide-context-specific neutral mutation rates are derived via the deconvolution of mutational signatures, performed by [deconstructSigs](https://github.com/raerose01/deconstructSigs). Next, all whole-exome and whole-genome tumor samples are used to estimate per-gene netural mutation rates using [dNdScv](https://github.com/im3sanger/dndscv), with support for tissue-specific covariates if supplied. Combining this information, the rate of neutral mutation at a particular SNV site in a tumor is that tumors's trinucleotide-context-specific mutation rate normalized by the neutral mutation rate of the gene. CancereffectsizeR uses maximum-likelihood methods to estimate the selection intensity at SNV sites by comparing the frequency of variants in the data with their expected frequency in the absence of selection.

## Version Note
[Version 0.1.0](https://github.com/Townsend-Lab-Yale/cancereffectsizeR/releases/tag/0.1.0) of this package was developed by Cannataro, V. L., Gaffney, S. G., and Townsend, J. P., as decribed in our 2018 JNCI paper [_Effect sizes of somatic mutations in cancer_](https://doi.org/10.1093/jnci/djy168). A user guide for v0.1.0 is available [here](https://github.com/Townsend-Lab-Yale/cancereffectsizeR/blob/master/user_guide/cancereffectsizeR_user_guide.md).

This README is for the current version, which has a number of improvements.

# Installation
CancereffectsizeR installation requires R 3.5.0 or later and a recent version of the devtools package.

```R
# If you don't have devtools, install it (or re-install if your version is ancient)
install.packages("devtools")

# Install cancereffectsizeR and dependencies
devtools::install_github("Townsend-Lab-Yale/cancereffectsizeR")

```

# Usage
Input data should be MAF format (either a text file or a data frame with the same format). If you have chronological tumor progression state information (e.g., stages 1-4, or primary/metastatic, or pre-/post-treatment), then you can provide this information to the analysis in order to have baseline site mutation rates and selection intensities calculated per state rather than all together.

```R

  # Create CESAnalysis object and define the chronological tumor progression states
  analysis = CESAnalysis(genome = "hg19", progression_order = c("Primary", "Metastatic"))
  
  # Load in an MAF and give the name of the column that identifies tumor progression state
  analysis = load_maf(analysis, maf = "my_wes_lung_adenocarcinoma_data.maf", progression_col = "Primary_Met")
 
  
  # Prepare for effect size analysis by calculating mutational signatures for each tumor,
  # running dN/dS to estimate gene-level mutation rates (in absence of selection),
  # and annotating MAF with gene and trinucleotide information.
  analysis = calc_baseline_mutation_rates(analysis, covariate_file = "lung_pca", cores = 4)
  
  # Calculate selection intensities
  # If you have multiple computing cores and the parallel library (and are not using Windows),
  # you can parallelize the operation by providing however many cores you have available.
  analysis = ces_snv(analysis, cores = 4)

  # Create a filtered table of selection intensities for all recurrently mutated sites
  results = snv_results(analysis)
  results = results[tumors_with_variant > 1]
```


