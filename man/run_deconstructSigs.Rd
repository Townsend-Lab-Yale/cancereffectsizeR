% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_deconstructSigs.R
\name{run_deconstructSigs}
\alias{run_deconstructSigs}
\title{cancereffectsizeR wrapper for deconstructSigs}
\usage{
run_deconstructSigs(
  tumor_trinuc_counts,
  signatures_df,
  signatures_to_remove,
  tri.counts.method
)
}
\arguments{
\item{tumor_trinuc_counts}{one-row data.frame of trinuc variant counts (in deconstructSigs order) for one tumor}

\item{signatures_df}{data.frame of signatures (see COSMIC v3 signatures included with package for format)}

\item{signatures_to_remove}{names of signatures in signatures_df to keep out of deconstructSigs and assign zero weights}

\item{tri.counts.method}{exome/genome trinucleotide content normalization argument to pass to deconstructSigs (see its docs)}
}
\value{
a data.frame of signature weights
}
\description{
This function gets called internally by trinuc_mutation_rates() for each tumor in a CESAnalysis, accepting
a data.frame of mutation counts and returning a data.frame of signature weights.
}
\keyword{internal}
