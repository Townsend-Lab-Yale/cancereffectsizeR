# Internal VCF parser

Used by vcfs_to_maf_table()

## Usage

``` r
read_vcf(vcf, sample_id, vcf_name = sample_id)
```

## Arguments

- vcf:

  VCF filename or VCF-like data.table.

- sample_id:

  1-length sample identifier.

- vcf_name:

  1-length identifier used in some user messages.

## Value

MAF-like data.table
