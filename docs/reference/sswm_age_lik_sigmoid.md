# Generalized sigmoid continuous-covariate SSWM likelihood

Constructs a negative log-likelihood function for site-level selection
under the strong-selection, weak-mutation (SSWM) model, with selection
intensity modeled as a four-parameter generalized sigmoid function of a
continuous covariate:

## Usage

``` r
sswm_age_lik_sigmoid(
  rates_tumors_with,
  rates_tumors_without,
  covariate_data,
  selection_results = NULL
)
```

## Arguments

- rates_tumors_with:

  vector of site-specific mutation rates for all tumors with variant

- rates_tumors_without:

  vector of site-specific mutation rates for all eligible tumors without
  variant

- covariate_data:

  A `data.table` containing `Unique_Patient_Identifier` and
  `covariate_value`. The `covariate_value` column contains the
  sample-level continuous covariate and must be coercible to numeric.
  For the generalized sigmoid model, covariate values must be strictly
  positive.

## Value

A function of a four-element numeric parameter vector corresponding to
`C`, `L`, `s`, and `m`. The returned function evaluates the negative
log-likelihood for the variant under the generalized sigmoid
continuous-covariate SSWM model.

## Details

\$\$ \gamma(x) = \exp(C) + (\exp(L)-\exp(C))\frac{x^s}{x^s + m^s}. \$\$

The parameters \\\exp(C)\\ and \\\exp(L)\\ define the lower and upper
asymptotic selection intensities, \\m\\ specifies the covariate value at
which selection intensity is halfway between the two asymptotes, and
\\s\\ is the shape exponent controlling how sharply the curve
transitions around \\m\\.

This function is a likelihood-function factory intended for use with
[`ces_variant_sigmoid()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant_sigmoid.md).
For each variant, baseline mutation rates and sample information are
supplied automatically by the cancer effect size inference workflow.

## See also

[`ces_variant_sigmoid`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant_sigmoid.md),
[`sswm_age_lik`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/sswm_age_lik.md),
[`sswm_age_lik_logistic`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/sswm_age_lik_logistic.md)

Other continuous selection likelihoods:
[`sswm_age_lik()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/sswm_age_lik.md),
[`sswm_age_lik_logistic()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/sswm_age_lik_logistic.md)
