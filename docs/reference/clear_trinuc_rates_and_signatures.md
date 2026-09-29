# Clear mutational signature attributions and related mutation rate information

Removes all data calculated or supplied via trinuc_mutation_rates,
set_signature_weights, set_trinuc_rates, etc. This function can be used
if you want to re-run signature analysis with different sample groupings
or parameters.

## Usage

``` r
clear_trinuc_rates_and_signatures(cesa)
```

## Arguments

- cesa:

  CESAnalysis
