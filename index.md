# (placeholder; leave)

<h1 style="margin-bottom:0px"><span class="grad1">cancer</span><span class="grad2">effect</span><span class="grad3">sizeR</span></h1>
#### Quantify somatic evolution in cancer

---

Welcome to cancereffectsizeR! This R package contains a variety of tools for analyzing somatic variant data and ultimately characterizing the evolutionary trajectories of cancers:

- Annotate somatic variants with built-in reference data, or create a [custom reference data set](articles/custom_refset_instructions.html) for almost any species/genome.
- Attribute mutations to mutational processes by extracting COSMIC signatures, or any custom set of signature definitions.
- Use provided tissue-specific covariates to inform calculation of gene mutation rates via [dNdScv](https://github.com/im3sanger/dndscv), or build your own [custom covariates](articles/create_custom_covariates.html).
- Test single-variant, stage/grade-specific, and epistatic models of selection, or define and test your own models.
- Arbitrarily batch variants by position, gene, or functional annotation, and quantify selection by batch.

We plan to extend cancereffectsizeR's functionality as we continue development, and we welcome your feedback, ideas, and bug reports.

## Installation and tutorial
For installation instructions and an overview of package features, see the [tutorial](articles/cancereffectsizeR.html). The quickstart section offers a condensed introduction to get you running a basic analysis in minutes.

## Selected publications
* **[Attribution of Cancer Origins to Endogenous, Exogenous, and Preventable Mutational Processes](https://academic.oup.com/mbe/article/39/5/msac084/6570859)**, *Molecular Biology and Evolution * (2022).<br>Variant effect estimates from cancereffectsizeR are employed in a novel method to determine the relative contributions of various mutational processes to oncogenesis.

* **[Effect Sizes of Somatic Mutations in Cancer](https://doi.org/10.1093/jnci/djy168)**, *Journal of the National Cancer Institute *(2018).<br>A pan-cancer analysis of cancer effects employs [Version 0.1.0](https://github.com/Townsend-Lab-Yale/cancereffectsizeR/releases/tag/0.1.0) of this package. This original version was developed by Vincent Cannataro, Stephen Gaffney, and Jeffrey Townsend.








