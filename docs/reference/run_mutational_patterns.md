# cancereffectsizeR wrapper for fit_to_signatures

This function gets called internally by trinuc_mutation_rates() for each
tumor in a CESAnalysis, accepting a data.frame of mutation counts and
returning fit_to_signatures output. Note: this function supports indels
if passed in the same format.

## Usage

``` r
run_mutational_patterns(
  tumor_trinuc_counts,
  signatures_df,
  signatures_to_remove,
  mp_strict_args = list(),
  bootstrap_mutations = FALSE
)
```

## Arguments

- tumor_trinuc_counts:

  matrix of trinuc variant counts where columns respond to tumors and
  order of trinucleotide changes match signatures_df

- signatures_df:

  data.frame of signatures (see COSMIC v3 signatures included with
  package for format)

- signatures_to_remove:

  names of signatures in signatures_df to keep out of MutationalPatterns
  and assign zero weights. Only occurs when strict == FALSE

- mp_strict_args:

  named list of additional arguments to fit_to_signatures_strict

- bootstrap_mutations:

  T/F (default FALSE) whether to run fit_to_signatures_bootstrapped()
  with n_boot = 1 instead of fit_to_signatures_strict().

## Value

a data.frame of signature weights
