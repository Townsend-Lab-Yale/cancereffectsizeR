# Sort and validate input variant IDs

Sorts input variant IDs by type, completes IDs by adding protein ID to
plain variant names (e.g. "KRAS G12C"), and ensures that IDs are valid
even if not present in annotations. This includes verifying that
reference alleles are correct in SNV IDs and that amino acid changes are
possible given the coding sequence.

## Usage

``` r
sort_and_validate_variant_ids(cesa, input_ids, drop_unannotated = FALSE)
```

## Arguments

- cesa:

  CESAnalysis

- input_ids:

  Variant names/IDs, typically from user.

- drop_unannotated:

  Whether to include variants that are not annotated in output.

## Value

List of variant_ids, with each element corresponding to one
variant_type.
