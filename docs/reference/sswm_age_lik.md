# Linear continuous-covariate SSWM likelihood

Constructs a negative log-likelihood function for site-level selection
under the strong-selection, weak-mutation (SSWM) model, with selection
intensity modeled as a linear function of a continuous covariate:

## Usage

``` r
sswm_age_lik(
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

## Value

A function of a numeric parameter vector `c(beta0, beta1)` that returns
the negative log-likelihood for the variant under the linear
continuous-covariate SSWM model.

## Details

\$\$\gamma(x) = \beta_0 + \beta_1 x,\$\$

where \\x\\ is the sample-level continuous covariate, \\\beta_0\\ is the
intercept, and \\\beta_1\\ is the covariate-dependent slope.

This function is a likelihood-function factory intended for use with
[`ces_variant_linear()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant_linear.md).
For each variant, the mutation rates and sample information required by
the likelihood are supplied automatically by the cancer effect size
inference workflow.

The returned function accepts a parameter vector `c(beta0, beta1)` and
returns the negative log-likelihood.

## See also

[`ces_variant_linear`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant_linear.md),
[`sswm_age_lik_logistic`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/sswm_age_lik_logistic.md),
[`sswm_age_lik_sigmoid`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/sswm_age_lik_sigmoid.md)

Other continuous selection likelihoods:
[`sswm_age_lik_logistic()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/sswm_age_lik_logistic.md),
[`sswm_age_lik_sigmoid()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/sswm_age_lik_sigmoid.md)
