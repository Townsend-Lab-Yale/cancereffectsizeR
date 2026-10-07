# Get SNVs that cause an amino acid change

An internal function to figure out the SNVs that can cause a given amino
acid substitution in a transcript

## Usage

``` r
aac_to_snv_ids(refcds_entry_name, aa_pos, aa_alt, bsg, refcds)
```

## Arguments

- aa_pos:

  Integer position of substitution on the transcript.

- aa_alt:

  Identity of substitution, either a three-letter code ("Lys") or "STOP"

- bsg:

  A BSgenome object for the genome build associated with the RefCDS
  entry

- refcds_entry:

  A RefCDS entry for the relevant transcript
