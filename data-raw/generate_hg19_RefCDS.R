
# Read in Gencode GTF (release 35, GRCh37-lifted version) a data table
gen = as.data.table(rtracklayer::import("~/Documents/reference/gencode.v35lift37.annotation.gtf"))

# Restrict to consensus coding sequences
ccds = gen[tag == "CCDS"]

# Hand-pick the two most significant CDKN2A transcripts and give them new gene names; discard other CDKN2A transcripts
# It's not generally recommended to rename genes; making an exception for historical reaons
ccds[protein_id == "ENSP00000307101.5", gene_name := "CDKN2A.p16INK4a"]
ccds[protein_id == "ENSP00000462950.1", gene_name :=  "CDKN2A.p14arf"]
ccds = ccds[gene_name != "CDKN2A"]


# One other tweak: specify an additional splice site position in the TP53 transcript that will be used
tp53_extra_splice_pos = list(ENSP00000269305.4 = 7579312)

output = build_RefCDS(gtf = ccds, genome = "hg19", additional_essential_splice_pos = tp53_extra_splice_pos)


refcds = output[[1]] 
gr_genes = output[[2]]

saveRDS(refcds, "RefCDS.rds")
saveRDS(gr_genes, "gr_genes.rds")
saveRDS(unique(gr_genes$names), "gene_names.rds")

