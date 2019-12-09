# takes at least several minutes to get all data from Ensembl
library(biomaRt)
mart = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", GRCh = 37)
attr = c("ensembl_gene_id", "ensembl_peptide_id", "chromosome_name", "strand", 
		"external_gene_name", "cds_length", "genomic_coding_start", "genomic_coding_end", 
		"cds_start", "cds_end")
cds_data = getBM(mart = mart, attributes = attr)
cds_data[cds_data == ""] = NA # a few fields return

# save data for personal use if frequently building hg19
# write.table(cds_data, "ensembl_hg19.txt", sep="\t", quote=F, row.names = F)

# pick a 1,000 "random" genes 
# need to get all transcripts from each gene, so you can't just randomly sample rows
set.seed(3019)
all_genes = unique(cds_data$ensembl_gene_id)
chosen = all_genes[sample(1:length(all_genes), 1000)]
cds_small = cds_data[cds_data$ensembl_gene_id %in% chosen,]

out_path = paste0(system.file(package="cancereffectsizeR"), 
	"/tests/testthat/testdata/ensembl_cds_hg19_small.rds")
saveRDS(cds_small, out_path)


## To generate the mini RefCDS object
genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
GenomeInfoDb::seqlevelsStyle(genome) = "NCBI" # script actually does this, anyway

# this will change soon (and need to handle TP53/CDKN2A, also)
refcds_and_gr_genes = ces_buildref(cds_small, genome)


