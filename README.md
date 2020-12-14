# <span style="color:#224570;font-size:115%"><em>cancereffectsizeR</em></span><br><span style="font-size:65%; color:#224570">Townsend Lab, Yale School of Public Health</span>

## Welcome
cancereffectsizeR aims to provide an intuitive and comprehensive set of tools to quantify the effects of somatic variants on cancer progression. Throughout the package, methods are designed to support as little or as much customization as you like:
- Annotate variants with built-in reference data, or create a custom reference data set for any genome build.
- Identify mutational processes in groups of samples by extracting COSMIC signatures, or any custom set of signature definitions.
- Use cancereffectsizeR's tissue-specific covariates to inform calculation of gene mutation rates (via [dNdScv](https://github.com/im3sanger/dndscv)), or supply your own covariates.
- Test single-variant, stage/grade-specific, and epistatic models of selection, or define and test your own models.
- Arbitrarily batch variants by position, gene, or functional annotation, and quantify selection by batch.

For a more detailed overview, see [Get Started](articles/cancereffectsizeR.html). Briefly, mutational signatures are extracted from each tumor's SNV mutation profile using [deconstructSigs](https://github.com/raerose01/deconstructSigs). The relative weights of biologically associated signatures are used to infer trinucleotide-context-specific relative rates of SNV mutations for each sample. Cohort-wide neutral gene mutation rates are calculated by calling [dNdScv](https://github.com/im3sanger/dndscv), ideally with tissue-specific covariates (included with cancereffectsizeR for many tissue types). Combining this information, the rate of neutral mutation at a particular variant site in a tumor is the context-specific mutation rate normalized by the gene mutation rate. Comparing rates of observed and expected mutation under a model of somatic selection allows an inference of selection intensity, which we also call cancer effect size.

Currently, cancereffectsizeR directly supports analysis of both noncoding and amino-acid-changing SNVs. Current methods provide some insight into interactions with other types of selection; for example, samples could be grouped by copy number status (or any other feature of interest) and assessed for differential SNV selection. We plan to extend cancereffectsizeR's functionality as we continue development, and we welcome your feedback, ideas, and bug reports.

## Installation
cancereffectsizeR requires R 3.5 or later and can be installed from its GitHub repository via the remotes package. Since we're not able to test all installation environments, please help by letting us know if you have any installation problems.

#### R 4.0 and later
On R 4.0 or later, run this:
```R
install.packages("remotes")
remotes::install_github("Townsend-Lab-Yale/cancereffectsizeR@*release")
```

#### R 3.5 and 3.6
Install as above, but first you need to install an archived version of a missing dependency:
```R
install.packages("https://cran.r-project.org/src/contrib/Archive/XML/XML_3.99-0.3.tar.gz", 
                 type = "source", repos = NULL)
```


## Publications and version note
Our 2018 JNCI paper [Effect sizes of somatic mutations in cancer](https://doi.org/10.1093/jnci/djy168) describes [Version 0.1.0](https://github.com/Townsend-Lab-Yale/cancereffectsizeR/releases/tag/0.1.0) of this package, which was developed by Cannataro, V. L., Gaffney, S. G., and Townsend, J. P. 

This site describes the latest version, which has methodological and performance improvements, including support for whole-genome and targeted gene sequencing data. A user guide for v0.1.0 is available [here](https://github.com/Townsend-Lab-Yale/cancereffectsizeR/blob/master/user_guide/cancereffectsizeR_user_guide.md).







