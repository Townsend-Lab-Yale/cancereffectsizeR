# step_selection_lik

Generates a log-likelihood function for step-specific selection: a
variant's scaled selection coefficient is allowed to differ across an
ordered sequence of tumor progression stages (for example, normal tissue
and primary tumor), rather than assuming one constant coefficient across
all samples (compare
[`sswm_lik()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/sswm_lik.md)).
All arguments to this likelihood function factory are automatically
supplied by
[`ces_variant_step()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant_step.md).

## Usage

``` r
step_selection_lik(
  rates_tumors_with,
  rates_tumors_without,
  sample_index,
  stage_mut_prop
)
```

## Arguments

- rates_tumors_with:

  vector of site-specific mutation rates for all tumors with variant

- rates_tumors_without:

  vector of site-specific mutation rates for all eligible tumors without
  variant

- sample_index:

  data.table with columns Unique_Patient_Identifier, group_index, and
  group_name, associating each sample with its (1-indexed) stage and the
  stage's display name. See
  [`assign_stage_index()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/assign_stage_index.md).

- stage_mut_prop:

  Numeric vector, in stage order, giving the proportion of the variant's
  gene-level mutation rate estimated to accumulate during each stage.
  Should sum to 1. See
  [`stage_mutation_proportions()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/stage_mutation_proportions.md).
