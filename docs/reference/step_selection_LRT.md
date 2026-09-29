# Likelihood ratio test for step-specific selection

Tests, for each variant (or compound variant), whether a step-specific
model of selection (see
[`ces_variant_step()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant_step.md))
fits significantly better than the constrained model in which one
selection intensity is shared across all stages.

## Usage

``` r
step_selection_LRT(
  cesa = NULL,
  step_run_name = NULL,
  simple_run_name = NULL,
  max_loglik = 1e+05
)
```

## Arguments

- cesa:

  CESAnalysis object containing the selection_results run(s).

- step_run_name:

  run_name of a previous
  [`ces_variant_step()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant_step.md)
  call.

- simple_run_name:

  Optional run_name of a previous
  [`ces_variant()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant.md)
  call over the same variants, to use instead of the built-in
  constrained fit (see details).

- max_loglik:

  Log-likelihoods above this value (default 1e5) are treated as
  implausible (usually a sign of optimizer non-convergence) and set to
  NA, along with the LRT statistic and p-value for that variant.

## Value

A data.table with columns variant_name, loglik_step, loglik_simple, df,
LRT_stat, and p_value.

## Details

By default, the constrained model is the one fit alongside the
step-specific model by
[`ces_variant_step()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant_step.md)
(`loglikelihood_constant` column). It uses the same likelihood and gene
mutation rates, so the two models are nested by construction.
Alternatively, supply `simple_run_name` to compare against a separate
[`ces_variant()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant.md)
run (with `return_fit = TRUE`) over the same variants. Only do this if
that run's model is nested in the step-specific model: each sample's
gene mutation rate in the
[`ces_variant()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant.md)
run must equal the cumulative rate for that sample's own stage.

The step-specific model has one additional free parameter per extra
stage beyond the first (that is, `num_stages - 1` more parameters than
the constrained model's single selection intensity), which sets the
test's degrees of freedom.
