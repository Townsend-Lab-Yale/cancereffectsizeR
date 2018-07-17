cancereffectsizeR user guide
================
Vincent Cannataro
2018-07-17

Introduction
============

`cancereffectsizeR` is an `R` package that may be used to calculate the effect size of single nucleotide variants (SNV) in cancer exome data[1]. It was designed for use with datasets in [MAF format](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) and works well with data from a large number of tumors, each with a large number of detected varaints. A previous version of `cancereffectsizeR` depended on output from `MutSigCV` (which runs in MATLAB) and can be found [here](https://github.com/Townsend-Lab-Yale/SNV_selection_intensity).

This vignette is broken into several sections detailing installation, usage, and a suite of `R` functions provided with this package that are useful for analyzing cancer exome data.

-   [Installing](#install)
-   [Calculating effect size](#effectsizecalc)
    -   [Preprocessing](#preprocess)
    -   [Calculating effect size](#effectsizecalc)

Installing
==========

`cancereffectsizeR` utilizes a few bioinformatic `R` packages, and thus requires the installation of some dependencies. If you already use `R` to bioinformatics, you likely already have these packages installed and you can directly install `cancereffectsizeR` in two lines.

``` r
install.packages("devtools",repos = "https://cloud.r-project.org")
devtools::install_github("Townsend-Lab-Yale/cancereffectsizeR")
```

However, if you just downloaded `R` today, you will need to install the dependencies.

``` r
source("https://bioconductor.org/biocLite.R")
biocLite("rtracklayer")
biocLite("GenomicRanges")
biocLite("BSgenome")
biocLite("BSgenome.Hsapiens.UCSC.hg19")
install.packages("deconstructSigs",repos = "https://cloud.r-project.org")
install.packages("devtools",repos = "https://cloud.r-project.org")
devtools::install_github("im3sanger/dndscv")
devtools::install_github("Townsend-Lab-Yale/cancereffectsizeR")
```

Calculating effect size
=======================

In this example, we will calculate the effect size of SNV within LGG (low grade glioma) using the TCGA LGG dataset provided at the [National Cancer Institute Genomic Data Common](https://portal.gdc.cancer.gov/). Specifically using the MAF generated with mutect variant calling (NCI UUID 2c0dab06-7af9-4380-bc83-13148298da19). Download and read the MAF into memory.

``` r
# MAF files are tab delim and contain 5 rows of header to skip
LGG_MAF <- read.delim(
  file = "../vignettes/TCGA.LGG.mutect.2c0dab06-7af9-4380-bc83-13148298da19.DR-7.0.somatic.maf",
  header = T,skip = 5,stringsAsFactors = F)
```

Preprocessing
-------------

### Converting hg38 to hg19

The latest TCGA data release has genomic variants in hg38 coordinates, so we need to convert these to hg19 so that the data is compatable with other packages utilized within `cancereffectsizeR`. `hg_converter` uses the [hg38ToHg19.over.chain](http://hgdownload.cse.ucsc.edu/gbdb/hg38/liftOver/hg38ToHg19.over.chain.gz) file and the `rtracklayer::liftOver` function to perform the conversion. Note that this function can convert between other builds with other \*.over.chain files.

``` r
# provide the path to the over.chain file. 
# I downloaded the file from 
# <http://hgdownload.cse.ucsc.edu/gbdb/hg38/liftOver/hg38ToHg19.over.chain.gz>
LGG_MAF <- hg_converter(chain = "~/Downloads/hg38ToHg19.over.chain",
                        maf_to_convert = LGG_MAF)
#> Loading in specified MAF...
#> Warning: closing unused connection 5 (~/Downloads/hg38ToHg19.over.chain)
#> Number of rows in the MAF that failed to convert:  4
```

### Adding columns with tumor sample identifier data and tumor allele data.

The TCGA does a great job documenting their data in a consistent fashion, so this step is optional if just using TCGA data. However, other sources may be less reliable, so these functions are provided in an attempt to get consistent tumor names and tumor allele nucleotides.

``` r
LGG_MAF <- unique_tumor_addition_function(MAF.file = LGG_MAF)
#> Summary statistics of the number of mutations per unique tumor:
#>     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#>     1.00    26.00    36.00    69.98    49.00 14335.00
```

``` r
LGG_MAF <- tumor_allele_adder(MAF = LGG_MAF)
```

### Removing DNV and TNV

We only want true single-nucleotide variant events in our data so we can get the best estimate of SNV mutation rate. Variant calling algorithms may mislabel di-nucleotide events as "SNV" and thus we need to remove these data before calculating mutation rate and effect size.

``` r
LGG_MAF <- DNP_TNP_remover(MAF = LGG_MAF)
#> Removing possible DNP
#> Total count of potential DNP removed:  168
#> DNP removal complete
```

Calculating cancer effect size
==============================

``` r
LGG_selection_output <- effect_size_SNV(MAF_file = LGG_MAF,
                                        covariate_file = "lgg_pca",
                                        genes_for_effect_size_analysis =
                                          c("IDH1","EGFR","TP53","BRAF"))
#> Checking if any reference alleles provided do not match reference genome...
#>     Note: 4 mutations removed for exceeding the limit of mutations per gene per sample
#> No mismatched mutations! Returning original input.
#> Removing all recurrent mutations...
#> Finding the number of mutations per tumor
#> Number of tumors over specified minimum mutation number of 50: 86
#> Cleaning input to only contain tumors above the minimum...
#> Calculating trinucleotide mutation counts...
#> Calculating individual tumor mutational signatures...
#> No id variables; using all as measure variables
#> Statistical summary of the proportion of the mutational signature in each tumor sample that is 'unknown'
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> 0.00000 0.02569 0.04938 0.06476 0.10100 0.19919
#> [1] Loading the environment...
#> [2] Annotating the mutations...
#>     Note: 1 samples excluded for exceeding the limit of mutations per sample
#>     Note: 4 mutations removed for exceeding the limit of mutations per gene per sample
#>     60 %...
#> [3] Estimating global rates...
#> [4] Running dNdSloc...
#> [5] Running dNdScv...
#>     Regression model for substitutions: all covariates were used (theta = 1.97).
#>     Regression model for indels: all covariates were used (theta = 0.848)
#> Calculating selection intensity...
```

``` r
print(head(
  LGG_selection_output$selection_output$all_mutations[
    order(
      LGG_selection_output$selection_output$all_mutations$selection_intensity,decreasing = T),],
  10),
  digits = 2)
#>     Gene AA_Pos Nucleotide_position Nuc_Ref Nuc_Change AA_Ref AA_Change
#> 4   IDH1    132                  NA    <NA>       <NA>      R         G
#> 1   IDH1    132                  NA    <NA>       <NA>      R         S
#> 3   IDH1    132                  NA    <NA>       <NA>      R         H
#> 9   EGFR    252                  NA    <NA>       <NA>      R         P
#> 66  TP53    175                  NA    <NA>       <NA>      R         G
#> 120 TP53    273                  NA    <NA>       <NA>      R         L
#> 2   IDH1    132                  NA    <NA>       <NA>      R         C
#> 93  TP53    230                  NA    <NA>       <NA>      T         P
#> 16  EGFR    324                  NA    <NA>       <NA>      R         L
#> 80  TP53    196                  NA    <NA>       <NA>      R         P
#>       gamma selection_intensity freq mutation_rate gene_AA_size dndscv_p
#> 4   5337133             2.2e+07   10       3.7e-09          415    0e+00
#> 1   1415599             5.8e+06    9       1.3e-08          415    0e+00
#> 3   3064357             3.6e+06  354       3.9e-07          415    0e+00
#> 9   2668101             2.8e+06    2       1.5e-09         1211    1e-14
#> 66   749607             1.3e+06    1       2.6e-09          394    0e+00
#> 120  491807             8.3e+05    4       1.6e-08          394    0e+00
#> 2    178880             7.1e+05   17       1.9e-07          415    0e+00
#> 93   314687             5.4e+05    1       6.3e-09          394    0e+00
#> 16   469748             5.0e+05    2       8.4e-09         1211    1e-14
#> 80   287301             4.9e+05    1       6.9e-09          394    0e+00
#>     dndscv_q Prop_tumors_with_specific_mut
#> 4    0.0e+00                        0.0197
#> 1    0.0e+00                        0.0177
#> 3    0.0e+00                        0.6969
#> 9    2.9e-11                        0.0039
#> 66   0.0e+00                        0.0020
#> 120  0.0e+00                        0.0079
#> 2    0.0e+00                        0.0335
#> 93   0.0e+00                        0.0020
#> 16   2.9e-11                        0.0039
#> 80   0.0e+00                        0.0020
#>     Prop_of_tumors_with_this_gene_mutated
#> 4                                   0.768
#> 1                                   0.768
#> 3                                   0.768
#> 9                                   0.061
#> 66                                  0.415
#> 120                                 0.415
#> 2                                   0.768
#> 93                                  0.415
#> 16                                  0.061
#> 80                                  0.415
```

[1] Cannataro, V. L., Gaffney, S. G., Townsend, J. P., “Effect sizes of somatic mutations in cancer" Preprint: <https://doi.org/10.1101/229724>
