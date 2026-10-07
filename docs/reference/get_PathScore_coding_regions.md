# Get PathScore coding regions

Returns GRanges that represent the coding sequence (CDS) definitions
used by PathScore. The hg19 version was created by running liftOver on
the hg38 intervals.

## Usage

``` r
get_PathScore_coding_regions(genome = "hg38")
```

## Arguments

- genome:

  Genome build: Either "hg39" (default) or "hg19".
