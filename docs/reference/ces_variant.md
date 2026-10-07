# Calculate cancer effects of variants

This function calculates variant effect sizes under the chosen model of
selection. Under the default model, a variant is assumed to have a
consistent scaled selection coefficient (cancer effect) across all
included samples.

## Usage

``` r
ces_variant(
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
  conf = 0.95
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

  Set to "basic" (default) or "sequential" (not yet available) to use
  built-in models of selection, or supply a custom function factory (see
  details).

- run_name:

  Optionally, a name to identify the current run.

- lik_args:

  Extra arguments, given as a list, to pass to custom likelihood
  functions.

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

  Cancer effect confidence interval width (NULL skips calculation,
  speeds runtime). Ignored when running custom models.

## Value

CESAnalysis object with selection results appended to the selection
output list

## Details

Definitions of the sample count columns in the effects output:

- included_with_variant: Number of samples that have the variant and
  were included in the inference.

- included_total: Number of samples that have coverage at the site and
  were included in the inference.

- held_out: Samples that have coverage at the site, but were held out of
  the inference due to `hold_out_same_gene_samples = TRUE`.

- uncovered: Samples that were not included in the inference because
  their sequencing did not cover the variant site.

Note that if a table of samples to include in the inference is specified
with `samples`, any CESAnalysis samples not present in the table will
not be included in any of the above accounts.

It's possible to pass in your own selection model. You'll need to create
a "function factory" that, for any variant, produces a likelihood
function that can be evaluated on the data. The first two arguments must
be `rates_tumors_with` and `rates_tumors_without`, which give the
baseline site mutation rates in samples with and without the variant.
The third argument must be `sample_index`, a data.table that associates
`Unique_Patient_Identifier` with group names and indices. (Your function
factory must accept this argument, but it doesn't have to use its
value.) Values for all three of these arguments will be calculated by
ces_variant and passed to your function factory automatically. Your
function can take whatever additional arguments you like, and you can
pass in values using `lik_args`. The likelihood function parameters that
ces_variant will optimize should be named and have default values. See
the source code of
[`sswm_lik()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/sswm_lik.md)
for an example.
