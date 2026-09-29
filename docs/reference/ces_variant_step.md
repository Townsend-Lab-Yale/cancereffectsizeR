# Calculate cancer effects of variants across ordered progression stages

Like
[`ces_variant()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant.md),
but under a step-specific model of selection: instead of assuming a
variant's scaled selection coefficient is constant across all included
samples, `ces_variant_step()` allows it to differ across an ordered
sequence of tumor progression stages (for example, normal tissue and
primary tumor, or an arbitrary number of stages).

## Usage

``` r
ces_variant_step(
  cesa = NULL,
  variants = select_variants(cesa, min_freq = 2),
  stage_mut_prop = NULL,
  stage_col = NULL,
  stage_order = NULL,
  sample_index = NULL,
  samples = character(),
  run_name = "auto",
  optimizer_args = list(method = "L-BFGS-B", lower = 0.001, upper = 1e+09),
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
  set to `[CESAnalysis]$variants`. Every variant (or compound variant)
  must map to exactly one gene present in `stage_mut_prop`.

- stage_mut_prop:

  A data.table of per-gene, per-stage mutation rate proportions, as
  produced by
  [`stage_mutation_proportions()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/stage_mutation_proportions.md).
  Its `p_<stage name>` columns must cover every stage named in
  `sample_index` (or implied by `stage_order`).

- stage_col:

  Name of a sample-level data column recording each sample's stage.
  Supply this and `stage_order` to have `sample_index` built
  automatically; alternatively, supply `sample_index` directly.

- stage_order:

  Values of `stage_col` in earliest-to-latest stage order (see
  [`assign_stage_index()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/assign_stage_index.md));
  used only when `sample_index` is not supplied directly.

- sample_index:

  A pre-built data.table associating samples with stages, as produced by
  [`assign_stage_index()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/assign_stage_index.md).
  Supply this, or `stage_col`/`stage_order`, but not both.

- samples:

  Which samples to include in inference. Defaults to all samples. Can be
  a vector of Unique_Patient_Identifiers, or a data.table containing
  rows from the CESAnalysis sample table. Every included sample must be
  covered by `sample_index`.

- run_name:

  Optionally, a name to identify the current run.

- optimizer_args:

  Named list of arguments to pass to the optimizer, bbmle::mle2.
  Defaults to bounding stage selection intensities to \[1e-3, 1e9\] via
  L-BFGS-B.

- return_fit:

  TRUE/FALSE (default FALSE): Embed model fit for each variant in a
  "fit" attribute of the selection results table. Use
  `attr(selection_table, 'fit')` to access the list of fitted models.
  Defaults to FALSE to save memory.

- hold_out_same_gene_samples:

  When finding likelihood of each variant, hold out samples that lack
  the variant but have any other mutations in the same gene. By default,
  TRUE when running with single variants, FALSE with a
  CompoundVariantSet.

- cores:

  Number of cores to use for processing variants in parallel (not useful
  for Windows systems).

- conf:

  Confidence interval width for stage selection intensities (NULL skips
  calculation, speeds runtime).

## Value

CESAnalysis object with selection results appended to the selection
output list

## Details

This function is a standalone counterpart to
[`ces_variant()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant.md)
(it does not call it, and
[`ces_variant()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant.md)
is unaffected by anything here): the underlying model of selection is
fundamentally different (multiple stage-specific selection coefficients,
rather than one), so the mechanics of confidence interval calculation
and output differ in ways not supported by
[`ces_variant()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant.md)'s
`model` argument.

Setting up a run requires two additional pieces of information beyond a
normal
[`ces_variant()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant.md)
call:

- Which stage each sample belongs to, in the form of a `sample_index`
  table (see
  [`assign_stage_index()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/assign_stage_index.md),
  or supply `stage_col`/`stage_order` here and it will be built for
  you).

- For each variant's gene, what proportion of the gene's mutation rate
  is estimated to accumulate during each stage, in the form of a
  `stage_mut_prop` table (see
  [`stage_mutation_proportions()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/stage_mutation_proportions.md)).

Because stage mutation-rate proportions are gene-level, each variant (or
compound variant) run through `ces_variant_step()` must map
unambiguously to exactly one gene.

The output selection table has one selection intensity column per stage,
named `si_<stage name>` (and, when `conf` is set, matching
`ci_low_<conf*100>_si_<stage name>`/`ci_high_<conf*100>_si_<stage name>`
columns), rather than the single `selection_intensity` column of
[`ces_variant()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant.md).
It does not include `included_with_variant`, `included_total`, or
`uncovered` columns, since sample accounting is inherently
stage-specific here.

Every included sample should carry the same gene mutation rates: the
final-stage cumulative rate, which the model divides among stages using
`stage_mut_prop` (see
[`stage_mutation_proportions()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/stage_mutation_proportions.md)).
A warning is issued if included samples belong to more than one gene
rate group.

Each variant is also fit under the constrained model in which all stages
share one selection intensity (`selection_intensity_constant` and
`loglikelihood_constant` columns). This constrained model is nested in
the step-specific model, using the same rates, so
[`step_selection_LRT()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/step_selection_LRT.md)
can test whether allowing selection to differ across stages
significantly improves fit. Use
[`plot_effects_step()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/plot_effects_step.md)
to visualize output.
