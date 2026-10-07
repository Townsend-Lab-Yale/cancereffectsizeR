# pairwise_epistasis_lik

For a pair of variants (or two groups of variants), creates a likelihood
function for a model of pairwise epistasis with a "strong mutation, weak
selection" assumption.

## Usage

``` r
pairwise_epistasis_lik(with_just_1, with_just_2, with_both, with_neither)
```

## Arguments

- with_just_1:

  two-item list of baseline rates in v1/v2 for tumors with mutation in
  just the first variant(s)

- with_just_2:

  two-item list of baseline rates in v1/v2 for tumors with mutation in
  just the second variant(s)

- with_both:

  two-item list of baseline rates in v1/v2 for tumors with mutation in
  both

- with_neither:

  two-item list of baseline rates in v1/v2 for tumors with mutation n
  neither

## Value

A likelihood function

## Details

The arguments to this function are automatically supplied by
[`ces_epistasis()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_epistasis.md)
and
[`ces_gene_epistasis()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_gene_epistasis.md).
