# Retrieve validated subset of CESAnalysis samples table

Retrieve validated subset of CESAnalysis samples table

## Usage

``` r
select_samples(cesa = NULL, samples = character())
```

## Arguments

- cesa:

  CESAnalysis

- samples:

  Vector of Unique_Patient_Identifiers, or data.table consisting of rows
  from a CESAnalysis samples table. If empty, returns full sample table.

## Value

data.table consisting of one or more rows from the CESAnalysis samples
table.
