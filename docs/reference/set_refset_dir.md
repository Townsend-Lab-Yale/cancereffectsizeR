# Set reference data directory

When working with custom reference data or loading a previously saved
CESAnalysis in a new environment, use this function to reassociate the
location of reference data with the analysis. (If
[`load_cesa()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/load_cesa.md)
didn't give you a warning when loading your analysis, you probably don't
need to use this function.)

## Usage

``` r
set_refset_dir(cesa, dir)
```

## Arguments

- cesa:

  CESAnalysis

- dir:

  path to data directory
