# clean_granges_for_cesa

Tries to format an input GRanges object to be compatible with a
CESAnalysis reference genome. Optionally, also applies padding to start
and end positions of ranges, stopping at chromosome ends. Either stops
with an error or returns a clean granges object.

## Usage

``` r
clean_granges_for_cesa(
  cesa = NULL,
  gr = NULL,
  padding = 0,
  refset_env = NULL,
  reduce_sort_strip = TRUE
)
```

## Arguments

- cesa:

  CESAnalysis; required unless refset_env is supplied directly.

- gr:

  GRanges object

- padding:

  How many bases to expand start and end of each position

- refset_env:

  A reference data set environment. Required if cesa is not specified.

- reduce_sort_strip:

  Unstrands and calls reduce(gr) which also drops all metadata columns,
  and then sorts the final gr. Default TRUE; use FALSE to preserve
  original gr structure.
