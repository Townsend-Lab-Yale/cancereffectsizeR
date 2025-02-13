% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trinuc_mutation_rates.R
\name{trinuc_snv_counts}
\alias{trinuc_snv_counts}
\title{Tabulate SNVs by trinucleotide context}
\usage{
trinuc_snv_counts(
  maf,
  genome,
  exclude_recurrent = FALSE,
  style = "MutationalPatterns"
)
}
\arguments{
\item{maf}{a cancereffectsizeR-style MAF data table}

\item{genome}{BSgenome reference genome (for looking up trinucleotide contexts)}

\item{exclude_recurrent}{Default FALSE. When TRUE, only mutations private to each sample are included in counts, in order to
reduce the influence of selection. (If you load more MAF data into the CESAnalysis later, recurrency may change.)}

\item{style}{"MutationalPatterns" or "deconstructSigs"}
}
\value{
Matrix or data frame of SNV counts, suitable for use with MutationalPatterns or
  deconstructSigs. Samples with zero passing SNVs will not appear.
}
\description{
This function produces trinucleotide-context-specific SNV counts from MAF data for
input to mutational signature extraction tools. Output can be tailored to meet
formatting requirements of MutationalPatterns or deconstructSigs, which are probably
similar to formats used by other tools.
}
