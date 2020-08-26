# <span style="color:#224570;font-size:115%">Welcome to <em>cancereffectsizeR</em></span><br><span style="font-size:65%; color:#224570">Townsend Lab, Yale School of Public Health</span>

cancereffectsizeR quantifies somatic selection for single-nucleotide variants in cancer cell lineages from tumor sequencing data. The package offers a straightforward, intuitive analytical framework that can be flexibly extended for customized analyses.

For an overview of core features, see [Get Started](articles/cancereffectsizeR.html). Briefly, trinucleotide-context-specific SNV mutation rates are calculated by extracting mutational signatures from each tumor's SNV mutation profile using [deconstructSigs](https://github.com/raerose01/deconstructSigs). Cohort-wide neutral gene mutation rates are calculated using [dNdScv](https://github.com/im3sanger/dndscv), with support for tissue-specific covariates. (Covariates data are included in cancereffectsizeR for many tissue types, and user-supplied data is also supported.) Combining this information, the rate of neutral mutation at a particular variant site in a tumor is the context-specific mutation rate normalized by the gene mutation rate. By comparing the observed frequency of variants with their expected frequencies, cancereffectsizeR infers selection intensities.

## Installation
Installing cancereffectsizeR requires R 3.5 or later and can be installed via the remotes or devtools packages. ()
```R
install.packages("devtools")
devtools::install_github("Townsend-Lab-Yale/cancereffectsizeR")
```

If you would like to contribute 



## Publications and version note
Our 2018 JNCI paper [Effect sizes of somatic mutations in cancer](https://doi.org/10.1093/jnci/djy168) describes [Version 0.1.0](https://github.com/Townsend-Lab-Yale/cancereffectsizeR/releases/tag/0.1.0) of this package, which was developed by Cannataro, V. L., Gaffney, S. G., and Townsend, J. P. 

This site describes the latest version, which has methodological and peformance improvements, including support for whole-genome and targeted gene sequencing data. A user guide for v0.1.0 is available [here](https://github.com/Townsend-Lab-Yale/cancereffectsizeR/blob/master/user_guide/cancereffectsizeR_user_guide.md).







