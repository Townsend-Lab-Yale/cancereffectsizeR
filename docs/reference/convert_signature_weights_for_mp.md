# Get MutationalPatterns contributions matrix

Reformat a signature weights table from mutational signature analysis
into the contributions matrix required for MutationalPatterns functions,
including visualizations.

## Usage

``` r
convert_signature_weights_for_mp(signature_weight_table)
```

## Arguments

- signature_weight_table:

  As created by trinuc_mutation_rates(); typically accessed via
  (CESAnalysis)\$mutational_signatures.
