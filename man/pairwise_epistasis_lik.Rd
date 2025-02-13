% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epistasis_objectives.R
\name{pairwise_epistasis_lik}
\alias{pairwise_epistasis_lik}
\title{pairwise_epistasis_lik}
\usage{
pairwise_epistasis_lik(with_just_1, with_just_2, with_both, with_neither)
}
\arguments{
\item{with_just_1}{two-item list of baseline rates in v1/v2 for tumors with mutation in just the first variant(s)}

\item{with_just_2}{two-item list of baseline rates in v1/v2 for tumors with mutation in just the second variant(s)}

\item{with_both}{two-item list of baseline rates in v1/v2 for tumors with mutation in both}

\item{with_neither}{two-item list of baseline rates in v1/v2 for tumors with mutation n neither}
}
\value{
A likelihood function
}
\description{
For a pair of variants (or two groups of variants), creates a likelihood function for a
model of pairwise epistasis with a "strong mutation, weak selection" assumption.
}
\details{
The arguments to this function are automatically supplied by \code{ces_epistasis()} and \code{ces_gene_epistasis()}.
}
