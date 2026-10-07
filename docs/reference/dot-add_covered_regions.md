# .add_covered_regions

.add_covered_regions

## Usage

``` r
.add_covered_regions(
  cesa,
  coverage_type,
  covered_regions,
  covered_regions_name,
  covered_regions_padding
)
```

## Arguments

- coverage_type:

  exome or targeted, if not using source_cesa

- covered_regions:

  A GRanges object or BED file path with genome build matching the
  target_cesa, if not using source_cesa

- covered_regions_name:

  A name to identify the covered regions, if not using source_cesa

- covered_regions_padding:

  optionally, add +/- this many bp to each interval in covered_regions

## Value

CESAnalysis given in target_cesa, with the new covered regions added
