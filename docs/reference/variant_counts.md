# Assess variant prevalence and coverage

Determine variant prevalence (and how many samples have sequencing
coverage) across your MAF data, or within different groups of samples.

## Usage

``` r
variant_counts(cesa, variant_ids = character(), by = character())
```

## Arguments

- cesa:

  CESAnalysis

- variant_ids:

  variant names ("KRAS G12C") or full variant IDs. If left empty, uses
  non-overlapping variants as returned by \`select_variants()\` with
  `min_freq = 1`.

- by:

  Optionally, a vector of one or more sample table columns. Variant
  prevalence and coverage data will be broken down by the groups defined
  by unique combinations of values in these columns.
