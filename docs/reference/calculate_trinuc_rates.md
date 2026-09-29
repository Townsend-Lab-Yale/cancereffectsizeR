# Calculate trinuc rates

Used internally to calculate trinuc rates from signature weights

## Usage

``` r
calculate_trinuc_rates(weights, signatures, tumor_names)
```

## Arguments

- weights:

  matrix of signature weights

- signatures:

  matrix of signatures

- tumor_names:

  names of tumors corresponding to rows of weights

## Value

matrix of trinuc rates where each row corresponds to a tumor

## Details

If any relative rate is less than 1e-9, we add the lowest
above-threshold rate to all rates and renormalize rates so that they sum
to 1.
