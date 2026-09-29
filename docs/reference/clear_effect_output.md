# Clear variant effect output

Remove output from previous ces_variant() runs from CESAnalysis

## Usage

``` r
clear_effect_output(cesa, run_names = names(cesa$selection))
```

## Arguments

- cesa:

  CESAnalysis

- run_names:

  Which previous runs to remove; defaults to removing all.
