# Estimate variant effects with a generalized sigmoid continuous-covariate model

Estimates cancer effects for individual or compound variants when the
selection intensity varies with a continuous sample-level covariate,
according to a four-parameter generalized sigmoid curve:

## Usage

``` r
ces_variant_sigmoid(
  cesa = NULL,
  variants = select_variants(cesa, min_freq = 2),
  samples = character(),
  model = "default",
  run_name = "auto",
  lik_args = list(),
  optimizer_args = if (identical(model, "default")) list(method = "L-BFGS-B", lower =
    0.001, upper = 1e+09) else list(),
  return_fit = FALSE,
  hold_out_same_gene_samples = "auto",
  cores = 1,
  conf = 0.95,
  optimizer = c("COBYLA", "ISRES"),
  constraint = 0.001
)
```

## Arguments

- cesa:

  CESAnalysis object

- variants:

  Which variants to estimate effects for, specified with a variant table
  such as from `[CESAnalysis]$variants` or
  [`select_variants()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/select_variants.md),
  or a `CompoundVariantSet` from
  [`define_compound_variants()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/define_compound_variants.md).
  Defaults to all recurrent mutations; that is,
  `[CESAnalysis]$variants[maf_prevalence > 1]`. To include all variants,
  set to `[CESAnalysis]$variants`.

- samples:

  Which samples to include in inference. Defaults to all samples. Can be
  a vector of Unique_Patient_Identifiers, or a data.table containing
  rows from the CESAnalysis sample table.

- model:

  A likelihood-function factory defining the selection model. The
  package provides likelihood functions for continuous-covariate
  selection that can be used directly; for the generalized sigmoid
  model, use
  [`sswm_age_lik_sigmoid()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/sswm_age_lik_sigmoid.md).
  A custom likelihood-function factory may also be supplied. See
  [`ces_variant()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant.md)
  for the general custom-model interface.

- run_name:

  Optionally, a name to identify the current run.

- lik_args:

  Named list of additional arguments passed to `model`. For
  continuous-covariate inference, this must include `covariate_data`, a
  table associating `Unique_Patient_Identifier` values with a numeric
  continuous covariate in `covariate_value`.

- optimizer_args:

  Named list of arguments to pass to the optimizer, bbmle::mle2. Use,
  for example, to choose optimization algorithm or parameter boundaries
  on custom models.

- return_fit:

  TRUE/FALSE (default FALSE): Embed model fit for each variant in a
  "fit" attribute of the selection results table. Use
  `attr(selection_table, 'fit')` to access the list of fitted models.
  Defaults to FALSE to save memory. Model fit objects can be of moderate
  or large size. If you run thousands of variants at once, you may
  exhaust your system memory.

- hold_out_same_gene_samples:

  When finding likelihood of each variant, hold out samples that lack
  the variant but have any other mutations in the same gene. By default,
  TRUE when running with single variants, FALSE with a
  CompoundVariantSet.

- cores:

  Number of cores to use for processing variants in parallel (not useful
  for Windows systems).

- conf:

  Nominal confidence level. Set to `NULL` to skip confidence interval
  calculation. For inference on continuous-model parameters, confidence
  intervals can also be calculated separately using the continuous-model
  confidence-interval routines.

- optimizer:

  Character scalar specifying the optimization algorithm. Currently
  supported values are `"COBYLA"` and `"ISRES"`, corresponding to the
  NLOPT algorithms `NLOPT_LN_COBYLA` and `NLOPT_GN_ISRES`, respectively.

- constraint:

  Positive numeric value specifying the minimum permitted separation
  between the two asymptotes. The optimization imposes
  `exp(L) - exp(C) >= constraint`.

## Value

A `CESAnalysis` object with a new entry appended to `selection_results`.
For each analyzed variant, the result contains the fitted generalized
sigmoid-model parameters ordered as `C`, `L`, `s`, and `m` and maximized
log-likelihood, together with variant annotations and applicable
sample-count information.

## Details

\$\$ \gamma(x) = \exp(C) + (\exp(L)-\exp(C))\frac{x^s}{x^s + m^s}. \$\$

The parameters \\\exp(C)\\ and \\\exp(L)\\ define the lower and upper
asymptotic selection intensities, \\m\\ specifies the covariate value at
which selection intensity is halfway between the two asymptotes, and
\\s\\ is the shape exponent controlling how sharply the curve
transitions around \\m\\.

Baseline mutation rates are obtained from the `CESAnalysis` object, and
model parameters are estimated separately for each variant by maximum
likelihood.

## See also

[`ces_variant`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant.md)

Other continuous selection models:
[`ces_variant_linear()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant_linear.md),
[`ces_variant_logistic()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant_logistic.md),
[`plot_effects_continuous()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/plot_effects_continuous.md)
