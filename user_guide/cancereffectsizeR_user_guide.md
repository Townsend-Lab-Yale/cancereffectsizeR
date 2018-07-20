cancereffectsizeR user guide
================
Vincent Cannataro
2018-07-20

-   [Introduction](#introduction)
-   [Installing](#installing)
-   [Calculating effect size](#calculating-effect-size)
    -   [Preprocessing](#preprocessing)
        -   [Converting hg38 to hg19](#converting-hg38-to-hg19)
        -   [Adding columns with tumor sample identifier data and tumor allele data.](#adding-columns-with-tumor-sample-identifier-data-and-tumor-allele-data.)
        -   [Removing DNV and TNV](#removing-dnv-and-tnv)
    -   [Calculating cancer effect size](#calculating-cancer-effect-size)
    -   [Interpreting results](#interpreting-results)
        -   [Selection intensity summary](#selection-intensity-summary)

Introduction
============

`cancereffectsizeR` is an `R` package that may be used to calculate the effect size of single nucleotide variants (SNV) in cancer exome data[1]. It was designed for use with datasets in [MAF format](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) and works well with data from a large number of tumors, each with a large number of detected varaints. A previous version of `cancereffectsizeR` depended on output from `MutSigCV` (which runs in MATLAB) and can be found [here](https://github.com/Townsend-Lab-Yale/SNV_selection_intensity).

This vignette is broken into several sections detailing installation, usage, and a suite of `R` functions provided with this package that are useful for analyzing cancer exome data.

<!-- * [Installing](#install)  -->
<!-- * [Calculating effect size](#effectsizecalc) -->
<!--     + [Preprocessing](#preprocess) -->
<!--     + [Calculating effect size](#effectsizecalc) -->
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
library(cancereffectsizeR) # load in the package 

# provide the path to the over.chain file. 
# I downloaded the file from 
# <http://hgdownload.cse.ucsc.edu/gbdb/hg38/liftOver/hg38ToHg19.over.chain.gz>
LGG_MAF <- hg_converter(chain = "~/Downloads/hg38ToHg19.over.chain",
                        maf_to_convert = LGG_MAF)
#> Warning in if (is.null(value)) {: closing unused connection 5 (~/Downloads/
#> hg38ToHg19.over.chain)
#> Loading in specified MAF...
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
------------------------------

The `effect_size_SNV()` function contains the entire pipeline necessary to calculate effect sizes, assuming your data is correctly preprocessed (in hg19 coordinates, no DNP, etc.). The function utilizes `deconstructSigs`[2] and `dndscv`[3], among other freely available `R` packages, to first determine the intrinsic mutation rate at all sites in each gene analyzed, and then the effect size given the detected prevalence of the mutation among sequenced tumors.

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

Interpreting results
--------------------

Output is grouped into four main sections.

-   Selection intensity / effect size information `LGG_selection_output$selection_output`
-   Gene-level mutation rate data `LGG_selection_output$mutation_rates`
-   Trinucleotide mutational signatures (from deconstructSigs) `LGG_selection_output$trinuc_data`
-   And, dndscv data (from dndscv) `LGG_selection_output$dndscvout`

### Selection intensity summary

The main output regarding effect size calculations is found within `LGG_selection_output$selection_output$all_mutations`. This dataframe contains useful information for each unique molecular variant in the analysis, included selection intensity, mutation rate, frequency of occurrence within the dataset (along with proportion of tumors with the specific variant (`Prop_tumors_with_specific_mut`) and proportion of tumors with a mutation in this gene (`Prop_of_tumors_with_this_gene_mutated`).

``` r
# selection intensity data

print(head(
  LGG_selection_output$selection_output$all_mutations[
    order(
      LGG_selection_output$selection_output$all_mutations$selection_intensity,decreasing = T),],
  6),
  digits = 2)
#>     Gene AA_Pos Nucleotide_position Nuc_Ref Nuc_Change AA_Ref AA_Change
#> 4   IDH1    132                  NA    <NA>       <NA>      R         G
#> 1   IDH1    132                  NA    <NA>       <NA>      R         S
#> 3   IDH1    132                  NA    <NA>       <NA>      R         H
#> 9   EGFR    252                  NA    <NA>       <NA>      R         P
#> 66  TP53    175                  NA    <NA>       <NA>      R         G
#> 120 TP53    273                  NA    <NA>       <NA>      R         L
#>       gamma selection_intensity prevalence_among_tumors mutation_rate
#> 4   5337133             2.2e+07                      10       3.7e-09
#> 1   1415599             5.8e+06                       9       1.3e-08
#> 3   3064357             3.6e+06                     354       3.9e-07
#> 9   2668101             2.8e+06                       2       1.5e-09
#> 66   749607             1.3e+06                       1       2.6e-09
#> 120  491807             8.3e+05                       4       1.6e-08
#>     gene_AA_size dndscv_p dndscv_q Prop_tumors_with_specific_mut
#> 4            415    0e+00  0.0e+00                        0.0197
#> 1            415    0e+00  0.0e+00                        0.0177
#> 3            415    0e+00  0.0e+00                        0.6969
#> 9           1211    1e-14  2.9e-11                        0.0039
#> 66           394    0e+00  0.0e+00                        0.0020
#> 120          394    0e+00  0.0e+00                        0.0079
#>     Prop_of_tumors_with_this_gene_mutated
#> 4                                   0.768
#> 1                                   0.768
#> 3                                   0.768
#> 9                                   0.061
#> 66                                  0.415
#> 120                                 0.415
```

A more detailed output is found within `LGG_selection_output$selection_output$complete_mutation_data`, which breaks down each molecular variants into the individual tumors it was found within.

``` r

print(head(LGG_selection_output$selection_output$complete_mutation_data),digits = 2)
#>   Gene Gene_size Percent_A Percent_T Percent_G Percent_C
#> 1 IDH1      1245      0.31      0.25      0.25       0.2
#> 2 IDH1      1245      0.31      0.25      0.25       0.2
#> 3 IDH1      1245      0.31      0.25      0.25       0.2
#> 4 IDH1      1245      0.31      0.25      0.25       0.2
#> 5 IDH1      1245      0.31      0.25      0.25       0.2
#> 6 IDH1      1245      0.31      0.25      0.25       0.2
#>   Nucleotide_Gene_Pos Nucleotide_chromosome_position Chromosome
#> 1                 395                      209113112          2
#> 2                 395                      209113112          2
#> 3                 395                      209113112          2
#> 4                 395                      209113112          2
#> 5                 395                      209113112          2
#> 6                 395                      209113112          2
#>   Reference_Nucleotide Alternative_Nucleotide Reference_Count
#> 1                    G                      A              NA
#> 2                    G                      A              NA
#> 3                    G                      A              NA
#> 4                    G                      A              NA
#> 5                    G                      A              NA
#> 6                    G                      A              NA
#>   Alternative_Count Tumor_origin Unique_patient_identifier
#> 1                NA TCGA-VM-A8CB              TCGA-VM-A8CB
#> 2                NA TCGA-DB-5273              TCGA-DB-5273
#> 3                NA TCGA-HT-7687              TCGA-HT-7687
#> 4                NA TCGA-DB-5279              TCGA-DB-5279
#> 5                NA TCGA-E1-5307              TCGA-E1-5307
#> 6                NA TCGA-QH-A6X4              TCGA-QH-A6X4
#>   Nucleotide_change_tally Nucleotide_mutation_rate Amino_acid_position
#> 1                     354                  3.9e-07                 132
#> 2                     354                  3.9e-07                 132
#> 3                     354                  3.9e-07                 132
#> 4                     354                  3.9e-07                 132
#> 5                     354                  3.9e-07                 132
#> 6                     354                  3.9e-07                 132
#>   Amino_acid_codon Codon_position Amino_acid_reference
#> 1              CGT              2                    R
#> 2              CGT              2                    R
#> 3              CGT              2                    R
#> 4              CGT              2                    R
#> 5              CGT              2                    R
#> 6              CGT              2                    R
#>   Amino_acid_alternative Amino_acid_mutation_rate Amino_acid_change_tally
#> 1                      H                  3.9e-07                     354
#> 2                      H                  3.9e-07                     354
#> 3                      H                  3.9e-07                     354
#> 4                      H                  3.9e-07                     354
#> 5                      H                  3.9e-07                     354
#> 6                      H                  3.9e-07                     354
#>     Gamma MAF_location Nucleotide_trinuc_context
#> 1 3064357           81                       CGT
#> 2 3064357          113                       CGT
#> 3 3064357          138                       CGT
#> 4 3064357          191                       CGT
#> 5 3064357          252                       CGT
#> 6 3064357          320                       CGT
#>   gene_level_synonymous_mutation_rate strand trinucs selection_intensity
#> 1                             1.6e-07      -    <NA>             3559278
#> 2                             1.6e-07      -    <NA>             3559278
#> 3                             1.6e-07      -    <NA>             3559278
#> 4                             1.6e-07      -    <NA>             3559278
#> 5                             1.6e-07      -    <NA>             3559278
#> 6                             1.6e-07      -    <NA>             3559278
```

`$selection_output` also contains information about the last gene analyzed, such as the nucleotide mutation rates at every position...

``` r
LGG_selection_output$selection_output$nucleotide_mutation_rates[,1:5]
#>   Pos. 1 Ref: A Pos. 2 Ref: T Pos. 3 Ref: G Pos. 4 Ref: G Pos. 5 Ref: C
#> A  0.000000e+00  2.384755e-08  9.587197e-08  1.206397e-07  7.138277e-09
#> T  2.470056e-08  0.000000e+00  4.841226e-08  3.493356e-08  2.266729e-07
#> G  4.467539e-08  1.520567e-08  0.000000e+00  0.000000e+00  1.256772e-09
#> C  1.022811e-08  7.883035e-08  1.557953e-08  1.380589e-08  0.000000e+00
```

... and the amino acid mutations at every position.

``` r
LGG_selection_output$selection_output$amino_acid_mutation_rates[,1:5]
#>      Pos. 1 Ref: Met Pos. 2 Ref: Ala Pos. 3 Ref: Ala Pos. 4 Ref: Leu
#> Phe     0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00
#> Leu     3.492867e-08    0.000000e+00    0.000000e+00    3.150253e-07
#> Ser     0.000000e+00    3.493356e-08    3.493356e-08    0.000000e+00
#> Tyr     0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00
#> Cys     0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00
#> Trp     0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00
#> Pro     0.000000e+00    1.380589e-08    1.380589e-08    5.587236e-08
#> His     0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00
#> Gln     0.000000e+00    0.000000e+00    0.000000e+00    3.289667e-08
#> Arg     1.520567e-08    0.000000e+00    0.000000e+00    1.991334e-08
#> Ile     1.598638e-07    0.000000e+00    0.000000e+00    0.000000e+00
#> Met     0.000000e+00    0.000000e+00    0.000000e+00    3.604021e-08
#> Thr     7.883035e-08    1.206397e-07    1.206397e-07    0.000000e+00
#> Asn     0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00
#> Lys     2.384755e-08    0.000000e+00    0.000000e+00    0.000000e+00
#> Val     4.467539e-08    2.266729e-07    2.266729e-07    1.342409e-08
#> Ala     0.000000e+00    2.127574e-07    2.350680e-07    0.000000e+00
#> Asp     0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00
#> Glu     0.000000e+00    7.138277e-09    7.138277e-09    0.000000e+00
#> Gly     0.000000e+00    1.256772e-09    1.256772e-09    0.000000e+00
#> STOP    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00
#>      Pos. 5 Ref: Ser
#> Phe     0.000000e+00
#> Leu     0.000000e+00
#> Ser     2.266729e-07
#> Tyr     0.000000e+00
#> Cys     2.952168e-08
#> Trp     0.000000e+00
#> Pro     0.000000e+00
#> His     0.000000e+00
#> Gln     0.000000e+00
#> Arg     2.700924e-08
#> Ile     3.604021e-08
#> Met     0.000000e+00
#> Thr     1.342409e-08
#> Asn     9.457594e-08
#> Lys     0.000000e+00
#> Val     0.000000e+00
#> Ala     0.000000e+00
#> Asp     0.000000e+00
#> Glu     0.000000e+00
#> Gly     4.801477e-08
#> STOP    0.000000e+00
```

[1] Cannataro, V. L., Gaffney, S. G., Townsend, J. P., “Effect sizes of somatic mutations in cancer" Preprint: <https://doi.org/10.1101/229724>

[2] Rosenthal, R., McGranahan, N., Herrero, J., Taylor, B. S., & Swanton, C. (2016). deconstructSigs: delineating mutational processes in single tumors distinguishes DNA repair deficiencies and patterns of carcinoma evolution. Genome Biology, 17(1), 31. <https://doi.org/10.1186/s13059-016-0893-4>

[3] Martincorena, I., Raine, K. M., Gerstung, M., Dawson, K. J., Haase, K., Van Loo, P., … Campbell, P. J. (2017). Universal Patterns of Selection in Cancer and Somatic Tissues. Cell. <https://doi.org/10.1016/j.cell.2017.09.042>
