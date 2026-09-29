# Variant-level pairwise epistasis

Calculate selection intensity under an assumption of pairwise epistasis
between pairs of variants. CompoundVariantSets are supported.

## Usage

``` r
ces_epistasis(
  cesa = NULL,
  variants = NULL,
  samples = character(),
  run_name = "auto",
  cores = 1,
  model = "default",
  lik_args = list(),
  pval_calc_fn = NULL,
  optimizer_args = list(),
  conf = 0.95,
  return_fit = FALSE
)
```

## Arguments

- cesa:

  CESAnalysis

- variants:

  To test pairs of variants, supply a list where each element is a
  2-length vector of CES-style variant IDs. Alternatively (and often
  more usefully), supply a CompoundVariantSet (see
  [`define_compound_variants()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/define_compound_variants.md))
  to test all pairs of compound variants in the set.

- samples:

  Which samples to include in inference. Defaults to all samples. Can be
  a vector of Unique_Patient_Identifiers, or a data.table containing
  rows from the CESAnalysis sample table.

- run_name:

  Optionally, a name to identify the current run.

- cores:

  number of cores for parallel processing of variant pairs

- model:

  Set to "default" to use built-in model of epistatic selection, or
  supply a custom function factory (see details).

- lik_args:

  Extra arguments, given as a list, to pass to custom likelihood
  functions.

- pval_calc_fn:

  For use with custom models; optional. A function that takes an
  epistasis model fit as input and returns p-values and other
  descriptives.

- optimizer_args:

  List of arguments to pass to the optimizer (bbmle::mle2).

- conf:

  confidence interval size from 0 to 1 (.95 -\> 95%); NULL skips
  calculation, reducing runtime.

- return_fit:

  TRUE/FALSE (default FALSE): Embed epistatic model fits for each
  variant pair in a "fit" attribute of the epistasis results table. Use
  `attr(my_results, 'fit')` to access the list of fitted models.

## Value

CESAnalysis with a table of epistatic inferences appended to list
`[CESAnalysis]$epistasis`. Some column definitions:

- variant_A, variant_B: Names for the two variants or merged sets of
  variants in each epistatic inference. For brevity in the case of
  merged variant sets, we say that a sample with any variant in variant
  set A "has variant A."

- ces_A0: Cancer effect (scaled selection coefficient) of variant A that
  acts in the absence of variant B.

- ces_B0: Cancer effect of variant B that acts in the absence of variant
  A.

- ces_A_on_B: Cancer effect of variant A that acts when a sample already
  has variant B.

- ces_B_on_A: Cancer effect of variant B that acts when a sample already
  has variant A.

- p_A_change: P-value of likelihood ratio test (LRT) that informs
  whether selection for variant A significantly changes after acquiring
  variant B. The LRT compares the likelihood of the full epistatic model
  to that of a reduced model in which ces_A0 and ces_A_on_B are set
  equal. The p-value is the probability, under the reduced model, of the
  likelihood ratio being greater than or equal to the ratio observed.

- p_B_change: P-value of likelihood ratio test (LRT) that informs
  whether selection for variant B significantly changes after acquiring
  variant A. The LRT compares the likelihood of the full epistatic model
  to that of a reduced model in which ces_B0 and ces_B_on_A are set
  equal. The p-value is the probability, under the reduced model, of the
  likelihood ratio being greater than or equal to the ratio observed.

- p_epistasis: P-value of likelihood ratio test that informs whether the
  epistatic model better explains the mutation data than a non-epistatic
  model in which selection for mutations in each variant are independent
  of the mutation status of the other variant. Quite often, p_epistasis
  will suggest a significant epistatic effect even though p_A_change and
  p_B_change do not suggest significant changes in selection for either
  variant individually. This is because the degree of co-occurrence can
  often be explained equally well by a strong change in selection for
  either variant.

- expected_nAB_epistasis: The expected number of samples with both A and
  B mutated under the fitted epistatic model. Typically, this will be
  very close to the actual number of AB samples (nAB).

- expected_nAB_null: The expected number of samples with both A and B
  mutated under a no-epistasis model.

- AB_epistatic_ratio: The ratio
  `expected_nAB_epistasis/expected_nAB_null`. Useful to gauge the
  overall impact of epistatic interactions on the co-occurrence of
  variants A and B. Since the expectations take mutation rates into
  account, this ratio is a better indicator than the relative
  frequencies of A0, B0, AB, 00 in the data set.

- nA0, nB0, nAB, n00: Number of (included) samples with mutations in
  just A, just B, both A and B, and neither.

- ces_A_null, ces_B_null: Cancer effects of A and B when estimated
  independently; that is, effects under a no-epistasis model.
