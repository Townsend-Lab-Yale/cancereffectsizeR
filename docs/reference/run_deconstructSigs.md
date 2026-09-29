# cancereffectsizeR wrapper for deconstructSigs

This function gets called internally by trinuc_mutation_rates() for each
tumor in a CESAnalysis, accepting a data.frame of mutation counts and
returning a data.frame of signature weights.

## Usage

``` r
run_deconstructSigs(
  tumor_trinuc_counts,
  signatures_df,
  signatures_to_remove,
  tri.counts.method
)
```

## Arguments

- tumor_trinuc_counts:

  one-row data.frame of trinuc variant counts (in deconstructSigs order)
  for one tumor

- signatures_df:

  data.frame of signatures (see COSMIC v3 signatures included with
  package for format)

- signatures_to_remove:

  names of signatures in signatures_df to keep out of deconstructSigs
  and assign zero weights

- tri.counts.method:

  exome/genome trinucleotide content normalization argument to pass to
  deconstructSigs (see its docs)

## Value

a data.frame of signature weights
