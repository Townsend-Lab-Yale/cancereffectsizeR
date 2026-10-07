# Tabulate SNVs by trinucleotide context

This function produces trinucleotide-context-specific SNV counts from
MAF data for input to mutational signature extraction tools. Output can
be tailored to meet formatting requirements of MutationalPatterns or
deconstructSigs, which are probably similar to formats used by other
tools.

## Usage

``` r
trinuc_snv_counts(
  maf,
  genome,
  exclude_recurrent = FALSE,
  style = "MutationalPatterns"
)
```

## Arguments

- maf:

  a cancereffectsizeR-style MAF data table

- genome:

  BSgenome reference genome (for looking up trinucleotide contexts)

- exclude_recurrent:

  Default FALSE. When TRUE, only mutations private to each sample are
  included in counts, in order to reduce the influence of selection. (If
  you load more MAF data into the CESAnalysis later, recurrency may
  change.)

- style:

  "MutationalPatterns" or "deconstructSigs"

## Value

Matrix or data frame of SNV counts, suitable for use with
MutationalPatterns or deconstructSigs. Samples with zero passing SNVs
will not appear.
