# Tissue-specific mutational signature exclusions

Get suggestions on signatures_to_remove for trinuc_mutation_rates for
COSMIC signatures v3 and later. For details, see [this
article](https://townsend-lab-yale.github.io/cancereffectsizeR/articles/cosmic_cancer_type_note.html)
on our website.

## Usage

``` r
suggest_cosmic_signature_exclusions(
  cancer_type = NULL,
  treatment_naive = NULL,
  quiet = FALSE
)
```

## Arguments

- cancer_type:

  See
  [here](https://townsend-lab-yale.github.io/cancereffectsizeR/articles/cosmic_cancer_type_note.html)
  for supported cancer type labels.

- treatment_naive:

  give TRUE if samples were taken pre-treatment; FALSE or leave NULL
  otherwise.

- quiet:

  (default false) for non-interactive use, suppress explanations and
  advice.

## Value

a vector of signatures to feed to the
[`trinuc_mutation_rates()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/trinuc_mutation_rates.md)
`signature_exclusions` argument.
