# add_covered_regions

add_covered_regions

## Usage

``` r
add_covered_regions(
  target_cesa = NULL,
  source_cesa = NULL,
  covered_regions = NULL,
  covered_regions_name = NULL,
  coverage_type = NULL,
  covered_regions_padding = 0
)
```

## Arguments

- target_cesa:

  CESAnalysis with annotated variants that the covered regions will be
  added to

- source_cesa:

  Another CESAnalysis to copy all covered regions from

- covered_regions:

  A GRanges object or BED file path with genome build matching the
  target_cesa, if not using source_cesa

- covered_regions_name:

  A name to identify the covered regions, if not using source_cesa

- coverage_type:

  exome, genome, or targeted (if not using source_cesa)

- covered_regions_padding:

  optionally, add +/- this many bp to each interval in covered_regions

## Value

CESAnalysis given in target_cesa, with the new covered regions added
