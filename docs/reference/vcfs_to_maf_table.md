# Read a VCF into an MAF-like table

This function loads VCF files into MAF-like tables. A
Tumor_Sample_Barcode column is added, and the contents of the
POS/REF/ALT fields are converted to match the style used by MAF files
for Start_Position/Reference_Allele/Tumor_Seq_Allele2. The VCF file
should represent high-confidence somatic variant calls for a single
tumor sample.

## Usage

``` r
vcfs_to_maf_table(vcfs, sample_ids)
```

## Arguments

- vcfs:

  Vector of VCF file paths, or a list of VCF-like data.tables, or a
  single data.table.

- sample_ids:

  Identifiers to populate Tumor_Sample_Barcode, one per VCF.

## Value

A single data.table with MAF-style fields, suitable for use with
cancereffectsizeR.
