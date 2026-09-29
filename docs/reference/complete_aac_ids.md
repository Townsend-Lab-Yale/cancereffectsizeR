# Create full AAC ID

For example, KRAS_G12C -\> KRAS_G12C_ENSP00000256078 (ces.refset.hg19).
In cases of multiple protein IDs per gene, will return more IDs than
input. Otherwise, input/output will maintain order.

## Usage

``` r
complete_aac_ids(partial_ids, refset)
```

## Arguments

- partial_ids:

  AAC variant id prefixes, such as "KRAS_G12C" or "MIB2 G395C"

- refset:

  reference data set (environment object)
