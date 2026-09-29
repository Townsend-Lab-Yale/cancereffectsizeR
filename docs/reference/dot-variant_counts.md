# Internal variant prevalence and coverage calculation

Called by variant_counts() (and select_variants()) with validated
inputs.

## Usage

``` r
.variant_counts(
  cesa,
  samples,
  snv_from_aac,
  noncoding_snv_id,
  by_cols = character()
)
```

## Arguments

- cesa:

  CESAnalysis

- samples:

  validated samples table

- snv_from_aac:

  data.table with columns aac_id, snv_id (validated and with annotations
  in CESAnalysis)

- noncoding_snv_id:

  vector of snv_ids to treat as noncoding variants

- by_cols:

  validated column names from sample table that are suitable to use for
  counting by.
