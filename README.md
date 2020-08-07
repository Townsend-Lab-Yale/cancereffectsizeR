# Welcome to cancereffectsizeR
cancereffectsizeR quantifies selection for somatic mutations in cancer from collections of tumor sequencing data. In brief, each tumor sample's relative rates of trinucleotide-context-specific substitution mutations are determined by inferring mutational processes via signature extraction with [deconstructSigs](https://github.com/raerose01/deconstructSigs). Next, all whole-exome and whole-genome samples are used to estimate per-gene neutral mutation rates using [dNdScv](https://github.com/im3sanger/dndscv), with support for tissue-specific covariates. Combining this information, the rate of neutral mutation at a particular variant site in a tumor is that tumor's trinucleotide-context-specific mutation rate normalized by the neutral mutation rate of the gene. cancereffectsizeR uses maximum-likelihood methods to estimate the selection intensity at substitution sites by comparing the frequency of variants in the data with their expected frequency in the absence of selection.

We are continuing to develop and improve cancereffectsizeR, so please don't hesitate to contact us with any questions or suggestions.

## Installation
Installing cancereffectsizeR requires R 3.5 or later and a recent version of the devtools package.
```R
install.packages("devtools")
devtools::install_github("Townsend-Lab-Yale/cancereffectsizeR")
```

## Publications and version note
Our 2018 JNCI paper [Effect sizes of somatic mutations in cancer](https://doi.org/10.1093/jnci/djy168) describes [Version 0.1.0](https://github.com/Townsend-Lab-Yale/cancereffectsizeR/releases/tag/0.1.0) of this package, which was developed by Cannataro, V. L., Gaffney, S. G., and Townsend, J. P. 

This site describes the latest version, which has methodological and peformance improvements, including support for whole-genome and targeted gene sequencing data. A user guide for v0.1.0 is available [here](https://github.com/Townsend-Lab-Yale/cancereffectsizeR/blob/master/user_guide/cancereffectsizeR_user_guide.md).







