# Introduction
CancereffectsizeR is an R package that estimates the selection intensity of somatic mutations in cancer from collections of tumor sequencing data. In brief, each tumor sample's trinucleotide-context-specific neutral mutation rates are derived via the deconvolution of mutational signatures, performed by [deconstructSigs](https://github.com/raerose01/deconstructSigs). Next, all whole-exome and whole-genome tumor samples are used to estimate per-gene netural mutation rates using [dNdScv](https://github.com/im3sanger/dndscv), with support for tissue-specific covariates if supplied. Combining this information, the rate of neutral mutation at a particular SNV site in a tumor is that tumors's trinucleotide-context-specific mutation rate normalized by the neutral mutation rate of the gene. CancereffectsizeR uses maximum-likelihood methods to estimate the selection intensity at SNV sites by comparing the frequency of variants in the data with their expected frequency in the absence of selection.

## Version Note
[Version 0.1.0](https://github.com/Townsend-Lab-Yale/cancereffectsizeR/releases/tag/0.1.0) of this package was developed by Cannataro, V. L., Gaffney, S. G., and Townsend, J. P., as decribed in our 2018 JNCI paper [Effect sizes of somatic mutations in cancer](https://doi.org/10.1093/jnci/djy168). A user guide for v0.1.0 is available [here](https://github.com/Townsend-Lab-Yale/cancereffectsizeR/blob/master/user_guide/cancereffectsizeR_user_guide.md).

This README is for the current version, which supports whole-genome and targeted gene sequencing data and has many improvements to performance and usability.

# Installation
CancereffectsizeR installation requires R 3.5.0 or later and a recent version of the devtools package.

```R
# If you don't have devtools, install it (or re-install if your version is ancient)
install.packages("devtools")

# Install cancereffectsizeR and dependencies
devtools::install_github("Townsend-Lab-Yale/cancereffectsizeR")

```

# Usage

## Create CESAanlysis
Begin by declaring your CESAnalysis object, the primary data structure of cancereffectsizeR. If you have chronological tumor progression state information (e.g., stages 1-4, primary/metastatic, or pre-/post-treatment), then you can provide this information to the analysis in order to have baseline site mutation rates and selection intensities calculated per state rather than all together. You will also supply your genome build.

```R
  # skip progression_order for a single-state selection model
  analysis = CESAnalysis(genome = "hg19", progression_order = c("Primary", "Metastatic")) 
```

Note: Currently only hg19 is supported, but detailed instructions are coming soon for generating custom reference data to run CES with any genome build. 

## Load MAF data
Next, load [MAF](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/#gdc-maf-format-v100) data from a text file or a data frame into a CESAnalysis object. If column names don't match MAF format specifications, you can supply your own column names. The only required columns are those specifying chromosome, position, reference and tumor alleles, and sample ID (Tumor_Sample_Barcode, by default); any other columns in your data will be ignored. If your CESAnalysis includes chronological tumor progression states (see Create CESAnalysis, above), also specify "progression_col". 

### Sequencing data coverage
By default, data is assumed to be derived from whole-exome sequencing. Whole-genome data and targeted sequencing data are also supported when the "coverage" option is specified. For targeted sequencing data, you must provided the set of target regions as a BED file or GRanges object. You can also supply this information for exome data, if you have it available, but it's not required. Run `?load_maf` for more information.

### LiftOver support
If the MAF data you are loading is from a different genome build than your CESAnalysis, you can use the "chain_file" option to supply a UCSC-style chain file, and then your MAF coordinates will be automatically converted with liftOver.

### Examples

```R
  ## Example 1
  analysis = CESAnalysis(genome = "hg19")
  
  # load in some whole-exome MAF data
  analysis = load_maf(analysis, maf = "wes_lung_adenocarcinoma.maf")
  
  # add some whole-genome data (and lift from hg38 to hg19 with a user-provided chain file)
  analysis = load_maf(analysis, maf = "wgs_luad.maf", coverage = "genome", chain_file = "hg38ToHg19.over.chain")
  
  # add some targeted gene sequencing data
  analysis = load_maf(analysis, maf = "cancer-gene-targeted_luad.maf", coverage = "targeted", 
                      covered_regions = "my_target_regions.bed", covered_regions_name = "TGS1")
  
  
  ## Example 2
  # Creating an analysis for data that has samples annotated as stage 1-4
  analysis = CESAnalysis(genome = "hg19", progression_order = c(1, 2, 3, 4))
  analysis = load_maf(analysis, maf = "my_multistage.maf", progression_col = "stage")
```


## Calculate mutation rates and estimate selection
After loading in data, run the following functions to prepare for calculation of cancer effect size. See the help messages for each function for advanced options.
- **trinucleotide_mutation_weights**: uses deconstructSigs to assign mutational signature weightings to each tumor sample, and from there calculates trinucleotide-context-specific SNV mutation rates
- **gene_level_mutation_rates**: uses dNdScv to calculate gene-level mutation rates. It's highly recommended to supply tissue-specific covariate data if available. For hg19, default covariates from dNdScv are also available.
- **annotate_gene_maf**: annotates MAF data with genome/transcriptome context and gene information

Finally, use **ces_snv** to find calculate effect sizes.

```R
  # Create CESAnalysis object and define the chronological tumor progression states
  analysis = CESAnalysis(genome = "hg19", progression_order = c("Primary", "Metastatic"))
  analysis = load_maf(analysis, maf = "luad_data.maf", progression_col = "pri_met")
  
  
  analysis = trinucleotide_mutation_weights(analysis)
  analysis = gene_level_mutation_rates(analysis, covariate_file = "lung_pca")
  analysis = annotate_gene_maf(analysis)
  
  
  # Calculate selection intensities
  analysis = ces_snv(analysis)

  # Create a filtered table of selection intensities for all recurrently mutated sites
  results = snv_results(analysis)
  results = results[tumors_with_variant > 1]
```

### Note on parallel processing
If you are running in a MacOS or Linux environment, you can speed up processing by providing multiple computing cores using the "cores" argument. You'll also need the R package parallel, which you can easily install with install.packages if you don't already have it. Run `parallel::detectCores()` to find out how many cores you have available.
```R
  # assuming you've already loaded MAF data into object "cesa"
  cesa = trinucleotide_mutation_weights(cesa, cores = 4)
  cesa = gene_level_mutation_rates(cesa, covariate_file = "breast_pca")
  cesa = annotate_gene_maf(cesa)
  cesa = ces_snv(cesa, cores = 4)
```



