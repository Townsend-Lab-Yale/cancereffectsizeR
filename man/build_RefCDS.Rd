% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_RefCDS.R
\name{build_RefCDS}
\alias{build_RefCDS}
\title{cancereffectsizeR's RefCDS builder}
\usage{
build_RefCDS(
  gtf,
  genome,
  use_all_transcripts = TRUE,
  cds_ranges_lack_stop_codons = TRUE,
  cores = 1,
  additional_essential_splice_pos = NULL,
  numcode = 1,
  chromosome_style = "NCBI"
)
}
\arguments{
\item{gtf}{Path of a Gencode-style GTF file, or an equivalently formatted data
table. See details for required columns (features). It's possible to build such a 
table using data pulled from biomaRt, but it's easier to use a GTF.}

\item{genome}{Genome assembly name (e.g., "hg19"); an associated BSgenome object must be
available to load. Alternatively, supply a BSgenome object directly.}

\item{use_all_transcripts}{T/F (default TRUE): Whether to use all complete transcripts or just the longest
one for each gene.}

\item{cds_ranges_lack_stop_codons}{The CDS records in Gencode GTFs don't include the
stop codons in their genomic intervals. If your input does include the stop 
codons within CDS records, set to FALSE.}

\item{cores}{how many cores to use for parallel computations}

\item{additional_essential_splice_pos}{Usually not needed. A list of
additional essential splice site positions to combine with those calculated
automatically by this function. Each element of the list should have a name
matching a protein_id in the input and consist of a numeric vector of
additional positions. This option exists so that mutations at chr17:7579312
on TP53 are treated as splice site mutations in cancereffectsizeR's default
hg19 reference data set. (Variants at this coding position, which are
always synonymous, have validated effects on splicing, even though the
position misses automatic "essential splice" annotation by 1 base.)}

\item{numcode}{(don't use) NCBI genetic code number; currently only code 1, the
standard genetic code, is supported}

\item{chromosome_style}{Chromosome naming style to use. Defaults to "NCBI". For the human
genome, that means 1, 2,..., as opposed to "UCSC" style (chr1, chr2, ...). Value gets
passed to genomeInfoDb's seqlevelsStyle().}
}
\value{
A two-item list: RefCDS (which is itself a big list, with each entry containing
  information on one coding sequence (CDS)), and a GRanges object that defines the
  genomic intervals covered by each CDS.
}
\description{
Based on the buildref function in Inigo Martincorena's package dNdScv, this function
takes in gene/transcript/CDS definitions and creates a dNdScv-style RefCDS object and
an associated GenomicRanges object also required to run dNdScv.
}
\details{
Required columns are seqnames, start, end, strand, gene_name, gene_id, protein_id, and
type.  Only rows that have type == "CDS" will be used. Strand should be
"+" or "-".

By default, only one the longest complete transcript is used from each gene in the
input. If you set use_all_transcripts = TRUE, then all complete transcripts will be
used, resulting in multiple RefCDS entries for some genes. If you do this, you may
want to first eliminate low-confidence or superfluous transcripts from the input data.
}
