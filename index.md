# (placeholder; leave)

<h1 style="margin-bottom:0px"><span class="grad1">cancer</span><span class="grad2">effect</span><span class="grad3">sizeR</span></h1>
#### Quantify somatic evolution in cancer

---

Welcome to cancereffectsizeR! This R package contains a variety of tools for analyzing somatic variant data and ultimately characterizing the evolutionary trajectories of cancers.

Features include:

- Annotate somatic variants with built-in reference data, or create a [custom reference data set](articles/custom_refset_instructions.html) for almost any species/genome.
- Attribute mutations to mutational processes by extracting COSMIC signatures, or any custom set of signature definitions.
- Use provided tissue-specific covariates to inform calculation of gene mutation rates via [dNdScv](https://github.com/im3sanger/dndscv), or build your own [custom covariates](articles/create_custom_covariates.html).
- Test single-variant, stage/grade-specific, and epistatic models of selection, or define and test your own models.
- Arbitrarily batch variants by position, gene, or functional annotation, and quantify selection by batch.


## Installation
cancereffectsizeR requires R 3.5 or later and can be installed from its GitHub repository via the remotes package. Since we're not able to test all installation environments, please help by letting us know if you have any installation problems. 

```R
options(timeout = 600)
install.packages("remotes")
remotes::install_github("Townsend-Lab-Yale/cancereffectsizeR@*release", 
                         dependencies = TRUE)
```
For a minimal installation, leave out the dependencies argument to install just required dependencies, rather than both required and suggested.

Regrettably, there is a bug in some older versions of the GenomeInfoDb package that may produce cryptic error messages in cancereffectsizeR, along the lines of `!anyNA(m32) is not TRUE`. If you encounter this issue, upgrade your Bioconductor version: 
```R
# Only necessary if the current Bioconductor version is <3.14.
if (BiocManager::version() < as.package_version('3.14')) {
  BiocManager::install(version = "3.14")`
}
```


## Methods overview
For a thorough introduction to the package's usage and underlying methods, see the [tutorial](articles/cancereffectsizeR.html). The quickstart section has a condesned intro to both usage and theory.

Currently, cancereffectsizeR directly supports analysis of both noncoding and amino-acid-changing SNVs. Current methods provide some insight into interactions with other types of selection; for example, samples could be grouped by copy number status (or any other feature of interest) and assessed for differential SNV selection. We plan to extend cancereffectsizeR's functionality as we continue development, and we welcome your feedback, ideas, and bug reports.



## Publications and version note
Our 2018 JNCI paper [Effect sizes of somatic mutations in cancer](https://doi.org/10.1093/jnci/djy168) describes [Version 0.1.0](https://github.com/Townsend-Lab-Yale/cancereffectsizeR/releases/tag/0.1.0) of this package, which was developed by Cannataro, V. L., Gaffney, S. G., and Townsend, J. P. 

This site describes the latest version, which has methodological and performance improvements, including support for whole-genome and targeted gene sequencing data. A user guide for v0.1.0 is available [here](https://github.com/Townsend-Lab-Yale/cancereffectsizeR/blob/main/user_guide/cancereffectsizeR_user_guide.md).







