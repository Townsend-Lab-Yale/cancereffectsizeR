# Assign samples to ordered progression stages

Builds the `sample_index` table required by
[`step_selection_lik()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/step_selection_lik.md)
(and, in turn,
[`ces_variant_step()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant_step.md)),
associating each sample with an integer index and a display name for its
position in an ordered sequence of tumor progression stages (for
example, normal tissue, primary tumor, metastasis).

## Usage

``` r
assign_stage_index(
  cesa = NULL,
  stage_col = NULL,
  stage_order = NULL,
  samples = character()
)
```

## Arguments

- cesa:

  CESAnalysis object

- stage_col:

  Name of a sample-level data column (as added via
  [`load_maf()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/load_maf.md)
  or
  [`load_sample_data()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/load_sample_data.md))
  that records each sample's stage.

- stage_order:

  Describes the values of `stage_col` in earliest-to-latest stage order.
  Either a character vector with one value of `stage_col` per stage
  (stage names default to these same values), or a named list where each
  element gives one-or-more values of `stage_col` that should be grouped
  into the same stage, with the element's name used as the stage's
  display name.

- samples:

  Which samples to include. Defaults to all samples in the CESAnalysis.
  Can be a vector of Unique_Patient_Identifiers, or a data.table
  containing rows from the CESAnalysis sample table.

## Value

A data.table with columns Unique_Patient_Identifier, group_index (1 =
earliest stage), and group_name.
