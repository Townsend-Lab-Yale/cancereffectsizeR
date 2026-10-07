# Save a CESAnalysis in progress

Saves a CESAnalysis to a file by calling using base R's saveRDS
function. Also updates run history for reproducibility. Files saved
should be reloaded with
[`load_cesa()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/load_cesa.md).

## Usage

``` r
save_cesa(cesa, file)
```

## Arguments

- cesa:

  CESAnalysis to save

- file:

  filename to save to (must end in .rds)

## Details

Note that the genome reference data associated with a CESAnalysis
(refset) is not actually part of the CESAnalysis, so it is not saved
here. (Saving this data with the analysis would make file sizes too
large.) When you reload the CESAnalysis, you can re-associate the
correct reference data.
