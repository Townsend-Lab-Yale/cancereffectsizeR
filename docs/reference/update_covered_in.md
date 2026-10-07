# update_covered_in

Updates the covered_in annotation for all variants to include all
covered regions in the CESAnalysis

## Usage

``` r
update_covered_in(cesa)
```

## Arguments

- cesa:

  CESAnalysis

## Value

CESAnalysis with regenerated covered-in annotations

## Details

Also updates internal cached output of select_variants().
