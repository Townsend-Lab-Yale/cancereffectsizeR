# sswm_lik

Generates log-likelihood function of site-level selection with "strong
selection, weak mutation" assumption. All arguments to this likelihood
function factory are automatically supplied by
[`ces_variant()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant.md).

## Usage

``` r
sswm_lik(rates_tumors_with, rates_tumors_without)
```

## Arguments

- rates_tumors_with:

  vector of site-specific mutation rates for all tumors with variant

- rates_tumors_without:

  vector of site-specific mutation rates for all eligible tumors without
  variant
