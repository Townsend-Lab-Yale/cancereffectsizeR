# (placeholder; leave)

<h1 style="margin-bottom:0px"><span class="grad1">cancer</span><span class="grad2">effect</span><span class="grad3">sizeR</span></h1>
#### Quantify somatic evolution in cancer

---

Welcome to cancereffectsizeR! This R package provides a variety of tools for analyzing somatic variant data and characterizing the evolutionary trajectories of cancers. Key package features and related theory are presented in a [recent article](https://aacrjournals.org/cancerres/article/83/4/500/716429/Estimation-of-Neutral-Mutation-Rates-and) in *Cancer Research*. As we continue development of cancereffectsizeR, we welcome your [feedback, questions, and bug reports](https://github.com/Townsend-Lab-Yale/cancereffectsizeR/issues).

## Installation, tutorial, and customizations
For simple installation instructions and a demonstration of a cancer effect analysis, see the [tutorial](articles/cancereffectsizeR.html). The quickstart section offers a condensed introduction to get you running a basic analysis in minutes.

The package has extensive support for testing specific research questions with customized analyses:

- Quantify the cancer effects of variants under the default model of selection, test epistatic models of selection, or define and test your own models. Arbitrarily batch variants by position, gene, or functional annotation, and quantify selection by batch. Assess differential selection among patient subgroups.
- Combining mutational signature analysis with cancer effect estimation, compare signatures' relative contributions to oncogenesis and mutagenesis, using either COSMIC signatures or any custom signature set.
- Annotate somatic variants with built-in reference data, or create a [custom reference data set](articles/custom_refset_instructions.html) for almost any species/genome.
- Use provided tissue-specific covariates to inform calculation of gene mutation rates via [dNdScv](https://github.com/im3sanger/dndscv), or build your own [custom covariates](articles/create_custom_covariates.html).

## Selected publications

* **[Estimation of Neutral Mutation Rates and Quantification of Somatic Variant Selection Using cancereffectsizeR](https://aacrjournals.org/cancerres/article/83/4/500/716429/Estimation-of-Neutral-Mutation-Rates-and)**, *Cancer Research* (2023).<br>This resource report discusses the package's key features and methods and presents analyses validating that cancer effects are a useful quantification of the cancer relevance of somatic variants.

* **[Attribution of Cancer Origins to Endogenous, Exogenous, and Preventable Mutational Processes](https://academic.oup.com/mbe/article/39/5/msac084/6570859)**, *Molecular Biology and Evolution*  (2022).<br>Cancer effects are incorporated into a novel method to calculate the relative contributions of various mutational processes to oncogenesis. Apply the method yourself with `mutational_signature_effects()`.

* **[Effect Sizes of Somatic Mutations in Cancer](https://doi.org/10.1093/jnci/djy168)**, *Journal of the National Cancer Institute *(2018).<br>A pan-cancer analysis of cancer effects employing [version 0.1.0](https://github.com/Townsend-Lab-Yale/cancereffectsizeR/releases/tag/0.1.0) of this package. This original version was developed by Vincent Cannataro, Stephen Gaffney, and Jeffrey Townsend.

## Citation
When reporting work that uses cancereffectsizeR, please cite our [resource report](https://aacrjournals.org/cancerres/article/83/4/500/716429/Estimation-of-Neutral-Mutation-Rates-and) published in *Cancer Research*:

>Mandell JD, Cannataro VL, Townsend JP. Estimation of neutral mutation rates and quantification of somatic variant selection using cancereffectsizeR. Cancer Research. 2023 Feb 15; 83(4):500-505. [doi:10.1158/0008-5472.CAN-22-1508](https://www.doi.org/10.1158/0008-5472.CAN-22-1508).









