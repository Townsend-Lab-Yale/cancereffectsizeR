## This GTF, or an updated version, will likely be used to generate version 2 of ces.refset.hg19.
## Version 1 had the same tweaks to CDKN2A and TP53, but the transcripts were pulled from biomaRt.
## Read in Gencode GTF (release 35, GRCh37-lifted version) as a data table

gen = as.data.table(rtracklayer::import("gencode.v35lift37.annotation.gtf"))

# Exome capture 
exome_bed = "xgen-exome-research-panel-targets.bed"

# Restrict to consensus coding sequences
ccds = gen[tag == "CCDS"]

# Hand-pick the two most significant CDKN2A transcripts and give them new gene names; discard other CDKN2A transcripts
# It's not generally recommended to rename genes; making an exception for historical reaons
ccds[protein_id == "ENSP00000307101.5", gene_name := "CDKN2A.p16INK4a"]
ccds[protein_id == "ENSP00000462950.1", gene_name :=  "CDKN2A.p14arf"]
ccds = ccds[gene_name != "CDKN2A"]


# One other tweak: specify an additional splice site position in the TP53 transcript that will be used
tp53_extra_splice_pos = list(ENSP00000269305.4 = 7579312)

refcds_output = build_RefCDS(gtf = ccds, genome = "hg19", additional_essential_splice_pos = tp53_extra_splice_pos)

# Call build_refset
create_refset(output_dir = "tmp",
              refcds_output = refcds_output,
              species_name = "human",
              genome_build_name = "hg19",
              BSgenome_name = "hg19", 
              supported_chr = c(1:22, 'X', 'Y'),
              default_exome_bed = exome_bed,
              exome_interval_padding = 100)
              


