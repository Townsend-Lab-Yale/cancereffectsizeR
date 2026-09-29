# get_ces_signature_set

For a given CES reference data collection and signature set name,
returns cancereffectsizeR's internal data for the signature set in a
three-item list: the signature set name, a data table of signature
metadata, and a signature definition data frame

## Usage

``` r
get_ces_signature_set(refset, name)
```

## Arguments

- refset:

  name of refset (if using a custom refset, it must be loaded into a
  CESAnalysis already)

- name:

  name of signature set
