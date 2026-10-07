# Assign pre-calculated regional mutation rates

This function allows you to specify regional rates of
mutation–calculated however you like–to samples in your CESAnalysis.
Rates can be assigned to all samples or to specified samples.

## Usage

``` r
set_gene_rates(
  cesa = NULL,
  rates = NULL,
  samples = character(),
  missing_genes_take_nearest = FALSE
)
```

## Arguments

- cesa:

  CESAnalysis object

- rates:

  A two-column data.table with either gene name or protein_id in column
  1 and rate in column 2

- samples:

  Which samples the input rates apply to. Defaults to all samples. Can
  be a vector of Unique_Patient_Identifiers, or a data.table containing
  rows from the CESAnalysis sample table.

- missing_genes_take_nearest:

  Set to TRUE to have each gene/protein_id missing from rates take the
  rate of the nearest non-missing gene/protein.

## Details

Provide rates in a data.table with two columns: gene name or protein ID
(character) and rate (numeric, non-negative). Gene names or protein IDs
must match those in CESAnalysis reference data. (Some reference data
sets, such as ces.refset.hg19, only allow gene-level rates.) If you
don't want to supply rates for every gene, set
`missing_genes_take_nearest = T` to have each missing gene or coding
region take the rate of the nearest non-missing one.
