# Load MAF somatic mutation data

Load MAF data from a text file or data table into your CESAnalysis.
Column names are expected to match MAF format specifications
(Chromosome, Start_Position, etc.). It's recommended to use
preload_maf() to prep the input (including, optionally, liftOver
conversion of genomic coordinates), but if you have clean MAF data, you
can run this function directly. By default, data is assumed to be
derived from whole-exome sequencing. Whole-genome data and targeted
sequencing data are also supported when the `coverage` option is
specified.

## Usage

``` r
load_maf(
  cesa = NULL,
  maf = NULL,
  maf_name = character(),
  coverage = "exome",
  covered_regions = NULL,
  covered_regions_name = NULL,
  covered_regions_padding = 0,
  sample_data_cols = character(),
  enforce_default_exome_coverage = FALSE
)
```

## Arguments

- cesa:

  CESAnalysis.

- maf:

  Path of tab-delimited text file in MAF format, or an MAF in data.table
  or data.frame format.

- maf_name:

  Optionally, a name to identify samples coming from the current MAF.
  Used to populate the maf_source field of the CESAnalysis samples
  table.

- coverage:

  exome, genome, or targeted (default exome).

- covered_regions:

  optional for exome, required for targeted: a GRanges object or a BED
  file of covered intervals matching the CESAnalysis genome.

- covered_regions_name:

  a name describing the covered regions (e.g.,
  "my_custom_targeted_regions"); required when covered_regions are
  supplied.

- covered_regions_padding:

  How many bases (default 0) to expand start and end of each
  covered_regions interval, to include variants called just outside of
  targeted regions. Consider setting from 0-100bp, or up to the
  sequencing read length. If the input data has been trimmed to the
  targeted regions, leave set to 0.

- sample_data_cols:

  MAF columns containing sample-level data (e.g., tumor grade) that you
  would like to have copied into the CESAnalysis samples table.

- enforce_default_exome_coverage:

  When loading default exome data, exclude records that aren't covered
  in the default exome capture intervals included with CES genome
  reference data (default FALSE).

## Value

CESAnalysis with the specified MAF data loaded. The MAF data table
includes CES-generated variant IDs, a list of all genes overlapping the
site, and top_gene and top_consequence columns that give the most
significant annotated coding changes for each mutation record.
Annotation precedence is determined by MAF prevalence (usually equal),
essential splice status, premature stop codon, nonsilent status, MAF
mutation prevalence across the transcript (often favors longer
transcripts), and finally alphabetical order. The columns are
recalculated when more data is loaded, so changes in MAF prevalence can
change which variants are highlighted. Note that
`[CESAnalysis]$variants` contains more information about all
top_consequence variants and all noncoding variants from the MAF.
