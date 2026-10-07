# Add sample data

Insert sample-level data into a CESAnalysis samples table.

## Usage

``` r
load_sample_data(cesa, sample_data)
```

## Arguments

- cesa:

  CESAnalysis

- sample_data:

  data.table or data.frame with a Unique_Patient_Identifier column to
  match the CESAnalysis samples table, with one row per sample. (It's
  okay if some samples aren't present in the table.)
