# Find samples with specified variants

A convenience function to identify samples with specific variants.

## Usage

``` r
samples_with(cesa, any_of = NULL)
```

## Arguments

- cesa:

  CESAnalysis

- any_of:

  Select samples with ANY of the given variant names/IDs, such as
  c("8:142506482_C\>G", "KRAS G12C"). When a gene has multiple
  transcripts in reference data, you may wish to use full IDs, such as
  "KRAS_G12C_ENSP00000256078".
