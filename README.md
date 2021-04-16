# <span style="color:#224570;font-size:115%"><em>cancereffectsizeR</em></span><br><span style="font-size:65%; color:#224570">Townsend Lab, Yale School of Public Health</span>

<!-- MarkdownTOC autolink="true" -->

- [Welcome](#welcome)
- [Getting started](#getting-started)
	- [Quick start, using Docker image](#quick-start-using-docker-image)
		- [Running container with Singularity](#running-container-with-singularity)
		- [Run package tests](#run-package-tests)
	- [Manual installation](#manual-installation)
		- [R 4.0 and later](#r-40-and-later)
		- [R 3.5 and 3.6](#r-35-and-36)
- [Publications and version note](#publications-and-version-note)

<!-- /MarkdownTOC -->

## Welcome
cancereffectsizeR aims to provide an intuitive and comprehensive set of tools to quantify the effects of somatic variants on cancer progression. Throughout the package, methods are designed to support as little or as much customization as you like:
- Annotate variants with built-in reference data, or create a [custom reference data set](articles/custom_refset_instructions.html) for any genome build.
- Identify mutational processes in groups of samples by extracting COSMIC signatures, or any custom set of signature definitions.
- Use cancereffectsizeR's tissue-specific covariates to inform calculation of gene mutation rates (via [dNdScv](https://github.com/im3sanger/dndscv)), or use your own [custom covariates](articles/create_custom_covariates.html).
- Test single-variant, stage/grade-specific, and epistatic models of selection, or define and test your own models.
- Arbitrarily batch variants by position, gene, or functional annotation, and quantify selection by batch.

For a more detailed overview, see [Get Started](articles/cancereffectsizeR.html). Briefly, mutational signatures are extracted from each tumor's SNV mutation profile using [deconstructSigs](https://github.com/raerose01/deconstructSigs). The relative weights of biologically associated signatures are used to infer trinucleotide-context-specific relative rates of SNV mutations for each sample. Cohort-wide neutral gene mutation rates are calculated by calling [dNdScv](https://github.com/im3sanger/dndscv), ideally with tissue-specific covariates (included with cancereffectsizeR for many tissue types). Combining this information, the rate of neutral mutation at a particular variant site in a tumor is the context-specific mutation rate normalized by the gene mutation rate. Comparing rates of observed and expected mutation under a model of somatic selection allows an inference of selection intensity, which we also call cancer effect size.

Currently, cancereffectsizeR directly supports analysis of both noncoding and amino-acid-changing SNVs. Current methods provide some insight into interactions with other types of selection; for example, samples could be grouped by copy number status (or any other feature of interest) and assessed for differential SNV selection. We plan to extend cancereffectsizeR's functionality as we continue development, and we welcome your feedback, ideas, and bug reports.


## Getting started

The sections below describe a [no-installation approach](#quick-start-using-docker-image) using Docker (or Singularity), and instructions on [how to install the package into your local R environment](#manual-installation).

### Quick start, using Docker image

The quickest way to try out cancereffectsizeR is to use the Docker container image on [Docker Hub](https://hub.docker.com/r/townsendlab/cancereffectsizer), which has the R package, dependencies, and R Studio bundled in. With Docker [installed and running](https://www.docker.com/get-started) on your machine, run the following command in the terminal to download the image and launch RStudio:
```shell script
docker run --rm -e PASSWORD=rstudiopassword -p 8787:8787 townsendlab/cancereffectsizer
```
RStudio will be accessible in your browser at http://localhost:8787, with username `rstudio` and custom password `rstudiopassword`.

To open RStudio with your current directory available as the default working directory (`/home/rstudio`), run:
```shell script
docker run --rm -e PASSWORD=rstudiopassword -v $PWD:/home/rstudio \
  -p 8787:8787 townsendlab/cancereffectsizer
```

Alternatively, to open the R interpreter, with your current directory as working directory, bound to `/work` in the container, run:
```shell script
docker run --rm -ti -v $PWD:/work -w /work townsendlab/cancereffectsizer R
```

#### Running container with Singularity

If you work in an HPC environment that uses [Singularity](https://sylabs.io/guides/latest/user-guide/) to run containers, you can launch RStudio with the following commands:
```shell script
export SINGULARITYENV_USER=$USER
export SINGULARITYENV_PASSWORD=rstudiopassword

mkdir -p /tmp/$USER run var-lib-rstudio-server
printf 'provider=sqlite\ndirectory=/var/lib/rstudio-server\n' > database.conf

singularity exec -C --home $(pwd):/home/rstudio \
  -B /tmp/${USER}:/tmp,$(pwd):/home/rstudio,run:/run \
  -B var-lib-rstudio-server:/var/lib/rstudio-server \
  -B database.conf:/etc/rstudio/database.conf \
  docker://townsendlab/cancereffectsizer rserver \
  --auth-none=0 --auth-pam-helper-path=pam-helper \
  --auth-timeout-minutes=0 --auth-stay-signed-in-days=30 \
  --www-address=0.0.0.0  --www-port 8787
```
As with running under Docker, RStudio will be accessible in your browser at http://localhost:8787, with username `rstudio` and custom password `rstudiopassword`. When running a remote machine, you will need to set up a tunnel for browser access, e.g.:
```shell script
ssh -N -L 8787:compute-node:8787 login-node
```

#### Run package tests

The container image includes the repository code at `/cancereffectsizeR`, which has a `tests` subdirectory. You can run the tests within the R interpreter or RStudio with:
```R
testthat::test_local('/cancereffectsizeR')
```


### Manual installation

cancereffectsizeR requires R 3.5 or later and can be installed from its GitHub repository via the remotes package. Since we're not able to test all installation environments, please help by letting us know if you have any installation problems.

#### R 4.0 and later
On R 4.0 or later, run this:
```R
install.packages("remotes")
options(timeout = 300)  # allow extra time to download BSgenome.Hsapiens.UCSC.hg19
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







