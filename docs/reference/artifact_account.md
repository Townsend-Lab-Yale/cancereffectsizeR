# Calculate relative rates of biological mutational processes

Sets artifact signature weights to zero and normalizes so that
biologically-associated weights sum to (1 - unattributed proportion) in
each sample.

## Usage

``` r
artifact_account(
  weights,
  signature_names,
  artifact_signatures = NULL,
  fail_if_zeroed = FALSE
)
```

## Arguments

- weights:

  data.table of signature weights (can have extra columns)

- signature_names:

  names of signatures in weights (i.e., all column names)

- artifact_signatures:

  vector of artifact signature names (or NULL)

- fail_if_zeroed:

  T/F on whether to exit if a tumor would have all-zero weights.
