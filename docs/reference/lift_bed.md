# Convert BED intervals between genome builds

Use this utility to convert BED intervals between genome coordinate
systems using liftOver. Only the chr/start/end fields of the input BED
are used (strand is ignored). The output GRanges will have no associated
seqinfo.

## Usage

``` r
lift_bed(bed, chain, outfile = NULL)
```

## Arguments

- bed:

  Pathname of a BED file, or a GRanges (typically loaded from a BED file
  with
  [`rtracklayer::import.bed()`](https://rdrr.io/pkg/rtracklayer/man/BEDFile-class.html)).

- chain:

  A UCSC-style chain file, or a Chain object (such as from
  [`rtracklayer::import.chain()`](https://rdrr.io/pkg/rtracklayer/man/Chain-class.html)).

- outfile:

  If not NULL, the returned GRanges will be saved to the specified path
  using `rtracklayer::export.chain()`.

## Value

GRanges representing lifted intervals from input `bed`.

## Details

A warning is given if the lifted intervals are less than 95% of the size
of the original intervals. When the BED input represents sequencing
target intervals, most of the input intervals will usually lift
successfully.
