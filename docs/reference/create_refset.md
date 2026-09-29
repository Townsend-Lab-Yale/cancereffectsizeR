# Create a custom refset

Use this function to create and save a directory of custom reference
data that can be used with cancereffectsizeR instead of supplied refsets
like `ces.refset.hg19`. All arguments are required except
default_exome/exome_interval_padding, which are recommended.

## Usage

``` r
create_refset(
  output_dir,
  refcds_dndscv,
  refcds_anno = NULL,
  species_name,
  genome_build_name,
  BSgenome_name,
  supported_chr = c(1:22, "X", "Y"),
  default_exome = NULL,
  exome_interval_padding = 0,
  transcripts = NULL,
  cores = 1
)
```

## Arguments

- output_dir:

  Name/path of an existing, writable output directory where all data
  will be saved. The name of this directory will serve as the name of
  the custom refset.

- refcds_dndscv:

  Transcript information in the two-item list (consisting of RefCDS and
  gr_genes) that is output by `build_RefCDS`. This transcript
  information will be used with dNdScv.

- refcds_anno:

  Transcript information in the two-item list (consisting of RefCDS and
  gr_genes) that is output by `build_RefCDS`. This transcript
  information will be used for cancereffectsizeR's annotations. If
  unspecified, the same reference information as supplied for dNdScv
  will be used.

- species_name:

  Name of the species, primarily for display (e.g., "human").

- genome_build_name:

  Name of the genome build, primarily for display (e.g., "hg19").

- BSgenome_name:

  The name of the BSgenome package to use (e.g., "hg19"); will used by
  cancereffectsizeR to load the reference genome via
  BSgenome::getBSgenome().

- supported_chr:

  Character vector of supported chromosomes. Note that cancereffectsizeR
  uses NCBI-style chromosome names, which means no chr prefixes ("X",
  not "chrX"). Mitochondrial contigs shouldn't be included since they
  would require special handling that hasn't been implemented.

- default_exome:

  A BED file or GRanges object that defines coding regions in the genome
  as might be used by an exome capture kit. This file (or GRanges) might
  be acquired or generated from exome capture kit documentation, or
  alternatively, coding regions defined in a GTF file (or the granges
  output by build_RefCDS()).

- exome_interval_padding:

  Number of bases to pad start/end of each covered interval, to allow
  for some variants to be called just outside of targeted regions, where
  there still may be pretty good sequencing coverage.

- transcripts:

  Additional information about coding (and, optionally, noncoding)
  transcripts from a Gencode GTF, supplied as a data.table. See the
  format provided in ces.refset.hg38. You'll have to match the format
  (including column names) pretty closely to get expected behavior.
  Noncoding transcripts are represented only by records with
  transcript_type = "transcript", and protein-coding transcripts are
  representing with transcript, CDS, and UTR records. Note that in
  Gencode format.

- cores:

  How many cores to use (default 1).

## Details

To run this function, you'll need to have output from
[`build_RefCDS()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/build_RefCDS.md).
