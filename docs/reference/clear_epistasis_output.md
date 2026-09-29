# Clear epistasis output

Remove previous epistatic effect estimations from CESAnalysis.

## Usage

``` r
clear_epistasis_output(cesa, run_names = names(cesa$epistasis))
```

## Arguments

- cesa:

  CESAnalysis.

- run_names:

  Which previous runs to remove; defaults to removing all.
