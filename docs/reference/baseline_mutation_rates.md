# Baseline mutation rate calculation

Calculates neutral mutation rates at specific sites based on gene
mutation rates and the relative trinucleotide-context-specific SNV
mutation rates of each sample

## Usage

``` r
baseline_mutation_rates(
  cesa,
  aac_ids = NULL,
  snv_ids = NULL,
  variant_ids = NULL,
  samples = character()
)
```

## Arguments

- cesa:

  CESAnalysis with gene mutation rates and tumor-specific
  trinucleotide-context-specific mutation rates already calculated

- aac_ids:

  vector of IDs for amino acid change variants

- snv_ids:

  vector of IDs for SNVs

- variant_ids:

  vector of mixed IDs (faster to use snv_ids and aac_ids for large jobs,
  if already known)

- samples:

  Which samples to calculate rates for. Defaults to all samples. Can be
  a vector of Unique_Patient_Identifiers, or a data.table containing
  rows from the CESAnalysis sample table.

## Value

a data table of mutation rates with one column per variant, and a
Unique_Patient_Identifier column identifying each row
