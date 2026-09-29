# Estimate selection under pairwise epistasis model

Estimate selection under pairwise epistasis model

## Usage

``` r
pairwise_variant_epistasis(
  cesa,
  variant_pair,
  samples,
  conf,
  compound_variants = NULL,
  model = "default",
  lik_args = list(),
  optimizer_args = list(),
  pval_calc_fn = NULL
)
```

## Arguments

- cesa:

  CESAnalysis

- variant_pair:

  2-length character of variant IDs, or 2-length numeric giving indices
  of CompoundVariantSet for the current two compound variants

- samples:

  Validated samples data.table (as from select_samples())

- compound_variants:

  If testing a pair of compound variants, the CompoundVariantSet
  defining them

- model:

  Passed from ces_epistasis or ces_gene_epistasis. Set to "default" to
  use built-in model of epistatic selection, or supply a custom function
  factory (see details).

- lik_args:

  Extra arguments, given as a list, passed from ces_epistasis or
  ces_gene_epistasis

- optimizer_args:

  List of arguments to pass to the optimizer (bbmle::mle2).

- pval_calc_fn:

  For use with custom models; optional. A function that takes an
  epistasis model fit as input and returns p-values and other
  descriptives.
