# Calculate per-gene, per-stage mutation rate proportions

For use with
[`step_selection_lik()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/step_selection_lik.md)/[`ces_variant_step()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant_step.md).
Takes cumulative per-gene mutation rates estimated separately for each
ordered progression stage (each from a
[`gene_mutation_rates()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/gene_mutation_rates.md)
call restricted to the samples at that stage; see the step-specific
selection vignette) and converts them into the proportion of a gene's
total mutation rate that is estimated to accumulate during each stage.

## Usage

``` r
stage_mutation_proportions(
  cesa = NULL,
  rate_cols = NULL,
  stage_names = NULL,
  rate_source = "dndscv_expected",
  on_invalid = "floor",
  floor_prop = 1e-06
)
```

## Arguments

- cesa:

  CESAnalysis object with gene rates already calculated once per stage,
  via repeated calls to
  [`gene_mutation_rates()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/gene_mutation_rates.md)
  (or
  [`set_gene_rates()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/set_gene_rates.md)),
  each restricted to the samples at that stage.

- rate_cols:

  Character vector of gene rate groups, one per stage, in
  earliest-to-latest order (for example,
  `c("rate_grp_1", "rate_grp_2")`). These name dNdScv runs in
  `cesa$dNdScv_results` (when `rate_source = "dndscv_expected"`) or
  columns of `get_gene_rates(cesa)` (when `rate_source = "gene_rates"`).
  Defaults to all `rate_grp_*` groups present, in ascending numeric
  order.

- stage_names:

  Optional character vector of display names for the stages, in the same
  order as `rate_cols`. Defaults to `rate_cols`.

- rate_source:

  `"dndscv_expected"` (default) for uncorrected rates calculated from
  full dNdScv output, or `"gene_rates"` to use `get_gene_rates(cesa)`
  as-is. See details.

- on_invalid:

  What to do when a gene's cumulative rates are not non-decreasing
  across stages: `"floor"` (default) floors negative/undefined stage
  contributions to `floor_prop` (of the gene's total rate) and rescales
  that gene's proportions to sum to 1, with a warning; `"error"` stops;
  `"NA"` sets all of that gene's proportions to `NA_real_`.

- floor_prop:

  Minimum stage proportion used when `on_invalid = "floor"` (default
  1e-6).

## Value

A data.table with a gene (or pid) identifier column, one proportion
column per stage (named `p_<stage_name>`), and one cumulative rate
column per stage (named `rate_<stage_name>`). Proportions in each row
sum to 1, except rows set to NA under `on_invalid = "NA"`.

## Details

By default (`rate_source = "dndscv_expected"`), each stage's cumulative
rate is dNdScv's covariate-based expected synonymous mutation count for
the gene (`exp_syn_cv`), divided by the gene's number of synonymous
sites and by the number of samples in the dNdScv run. These are
uncorrected rates: unlike the rates that
[`gene_mutation_rates()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/gene_mutation_rates.md)
assigns, they are not adjusted toward each gene's observed synonymous
mutation count. Using the same kind of estimate at every stage keeps
stage proportions internally consistent (the adjustment is skipped
altogether when dNdScv's overdispersion parameter theta is below 1,
which is common for sparse early-stage samples). This requires full
dNdScv output, so run each stage's
[`gene_mutation_rates()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/gene_mutation_rates.md)
call with `save_all_dndscv_output = TRUE`. Alternatively,
`rate_source = "gene_rates"` uses the rates in `get_gene_rates(cesa)`
as-is (for example, rates supplied with
[`set_gene_rates()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/set_gene_rates.md)).

[`ces_variant_step()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant_step.md)
expects every sample to carry the same, final-stage cumulative gene
rate, which the step-specific model then divides among stages using
these proportions. After calculating proportions, clear the per-stage
rates with
[`clear_gene_rates()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/clear_gene_rates.md)
and assign the final-stage cumulative rate (the last `rate_<stage name>`
column of the output) to all samples with
[`set_gene_rates()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/set_gene_rates.md).

Cumulative rates should be non-decreasing across stages, since later
stages by definition cover at least as much mutational "time" as earlier
ones. When dN/dS-based rate estimates are noisy, a later stage's
estimated cumulative rate can come out below an earlier one for a given
gene, which would produce an invalid (negative) proportion. When this
happens, affected genes are flagged with a warning, and handled
according to `on_invalid`.
