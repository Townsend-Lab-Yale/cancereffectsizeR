# Ensure AAC IDs are valid for a given reference data set

Given a vector of AAC IDs, determines whether each is valid/possible,
returning NULL if all are valid, or a list of problems, or an error if
input isn't parseable.

## Usage

``` r
validate_aac_ids(aac_ids, refset)
```

## Arguments

- aac_ids:

  AAC variant IDs

- refset:

  reference data set (environment object)

## Details

An ID is invalid if the gene and/or protein_id are not in the reference
data, or if the reference amino acid ("G" in KRAS_G12C) is incorrect, or
if there is no possible SNV that can create the proposed change. For
example, KRAS G12C is possible when the codon (GGT) acquires a G\>T
substitution in the first position, while G12K is not possible because
no single substitution can transform GGT into a lysine codon.
