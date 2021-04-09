# Partly based on https://github.com/Bioconductor/bioconductor_docker/blob/master/Dockerfile
FROM rocker/rstudio:4.0.5

# nuke cache dirs before installing pkgs; tip from Dirk E fixes broken img
RUN rm -f /var/lib/dpkg/available && rm -rf  /var/cache/apt/*

# issues with '/var/lib/dpkg/available' not found
# this will recreate
RUN dpkg --clear-avail

# This is to avoid the error 'debconf: unable to initialize frontend: Dialog'
ENV DEBIAN_FRONTEND noninteractive

ARG CRAN_URL="https://mran.microsoft.com/snapshot/2021-04-07"
ARG BIOC_URL="http://bioconductor.org/packages/3.12/bioc"
ARG N_CPU_BUILD=1

RUN apt-get update \
	&& apt-get install -y --no-install-recommends apt-utils \
	&& apt-get install -y --no-install-recommends \
	## Basic deps
	gdb \
	libxml2-dev \
	python3-pip \
	libz-dev \
	liblzma-dev \
	libbz2-dev \
	libpng-dev \
	libgit2-dev \
	&& rm -rf /var/lib/apt/lists/* \
	&& echo "options(repos = c(CRAN = '${CRAN_URL}', BIOC = '${BIOC_URL}'), download.file.method = 'libcurl')" \
		>> /usr/local/lib/R/etc/Rprofile.site \
	# cran, excl deconstructSigs
	&& install2.r -s -e -n ${N_CPU_BUILD} BH BiocManager MASS Matrix R6 RColorBrewer RCurl Rcpp XML ade4 askpass base64enc \
		bbmle bdsmatrix bitops brio cachem callr cli clipr colorspace cpp11 crayon credentials curl data.table \
		desc digest downlit ellipsis evaluate fansi farver fastmap formatR fs futile.logger futile.options gert \
		ggplot2 ggrepel gh gitcreds glue gtable highr hms htmltools httr ini isoband jsonlite knitr labeling \
		lambda.r lattice lazyeval lifecycle magrittr markdown matrixStats memoise mgcv mime mockr munsell \
		mvtnorm nlme numDeriv openssl pbapply pillar pixmap pkgconfig pkgdown plyr poilog prettyunits processx \
		progress ps purrr ragg rappdirs rematch2 renv reshape2 rlang rmarkdown rprojroot rstudioapi scales \
		segmented seqinr snow sp stringi stringr sys systemfonts textshaping tibble tinytex usethis utf8 vctrs \
		viridisLite whisker withr xfun xml2 yaml zip remotes \
	&& /usr/local/lib/R/site-library/littler/examples/installBioc.r BSgenome BSgenome.Hsapiens.UCSC.hg19 Biobase \
		BiocGenerics BiocParallel BiocVersion Biostrings DelayedArray GenomeInfoDb GenomeInfoDbData \
		GenomicAlignments GenomicRanges IRanges MatrixGenerics Rhtslib Rsamtools S4Vectors SummarizedExperiment \
		XVector rtracklayer zlibbioc deconstructSigs \
	&& installGithub.r -u FALSE im3sanger/dndscv "Townsend-Lab-Yale/cancereffectsizeR@*release" \
	&& installGithub.r -u FALSE "Townsend-Lab-Yale/ces.refset.hg19@*release"
