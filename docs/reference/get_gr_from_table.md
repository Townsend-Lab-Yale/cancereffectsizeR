# Get GRanges from chr/start/end table

Mainly built for select_variants() output, and uses the center_nt_pos on
AACs (rather than all from start-end). Assumes MAF-like coordinates
(1-based, closed).

## Usage

``` r
get_gr_from_table(variant_table)
```

## Arguments

- variant_table:

  data.table
