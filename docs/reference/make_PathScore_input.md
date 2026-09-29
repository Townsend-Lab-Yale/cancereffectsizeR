# Make a PathScore input file from MAF data

Produce a file from MAF data for use with PathScore, a web tool for
identifying significantly altered pathways in cancer. See
https://pathscore.publichealth.yale.edu/ and
https://doi.org/10.1093/bioinformatics/btw512.

## Usage

``` r
make_PathScore_input(maf, file = NULL, genome = "hg38")
```

## Arguments

- maf:

  An MAF-like data.table, as from preload_maf(). If the MAF has a column
  named `annot`, this column will be preserved in output (since
  PathScore supports an optional annotation column with this name).

- file:

  Name/path where PathScore-formatted data should be written. Will be a
  tab-delimited text file.

- genome:

  The genome build associated with the MAF file. Must be "hg38"
  (default) or "hg19".

## Value

A copy of the PathScore-formatted data as a data.table.
