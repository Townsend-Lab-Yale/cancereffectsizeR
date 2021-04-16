# Partly based on https://github.com/Bioconductor/bioconductor_docker/blob/master/Dockerfile
# Bioconductor 3.12 built for R 4.0.3
FROM rocker/rstudio:4.0.3

# nuke cache dirs before installing pkgs; tip from Dirk E fixes broken img
RUN rm -f /var/lib/dpkg/available && rm -rf  /var/cache/apt/*

# issues with '/var/lib/dpkg/available' not found
# this will recreate
RUN dpkg --clear-avail

# This is to avoid the error 'debconf: unable to initialize frontend: Dialog'
ENV DEBIAN_FRONTEND noninteractive

ARG N_CPU_BUILD=1

RUN apt-get update \
	&& apt-get install -y --no-install-recommends apt-utils \
	&& apt-get install -y --no-install-recommends \
	libxml2-dev \
	libz-dev \
	liblzma-dev \
	libbz2-dev \
	libgit2-dev \
	&& rm -rf /var/lib/apt/lists/* \
	&& install2.r -s -e -n ${N_CPU_BUILD} BiocManager remotes mockr devtools \
	&& echo "options(timeout = 300)" >> /usr/local/lib/R/etc/Rprofile.site \
	&& /usr/local/lib/R/site-library/littler/examples/installBioc.r \
		BSgenome.Hsapiens.UCSC.hg19 deconstructSigs \
	&& installGithub.r -u FALSE "im3sanger/dndscv"

COPY . /cancereffectsizeR

RUN R -e 'remotes::install_local("/cancereffectsizeR", upgrade = FALSE, dependencies = TRUE)'

RUN installGithub.r -u FALSE -d FALSE "Townsend-Lab-Yale/ces.refset.hg19@*release"
