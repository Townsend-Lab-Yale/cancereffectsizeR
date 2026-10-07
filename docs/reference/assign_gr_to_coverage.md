# assign_gr_to_coverage

Adds a validated GRanges object as a CESAnalysis's coverage set. Called
by add_covered_regions() after various checks pass.

## Usage

``` r
assign_gr_to_coverage(cesa, gr, covered_regions_name, coverage_type)
```

## Arguments

- cesa:

  CESAnalysis to receive the gr

- gr:

  GRanges

- covered_regions_name:

  unique name for the covered regions

- coverage_type:

  "exome" or "targeted"

## Details

Special handling occurs if covered_regions_name is "exome+".
