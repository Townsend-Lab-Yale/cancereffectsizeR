prev_dir = setwd(system.file("tests/test_data/", package = "cancereffectsizeR"))

## Previous method for getting CDS data using biomaRT
## The output ensembl_hg19.txt could still be used with the new build_RefCDS function
## with some pretty simple manipulations, but it's now recommend to just use a GTF instead
# library(biomaRt)
# mart = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", GRCh = 37)
# attr = c("ensembl_gene_id", "ensembl_peptide_id", "chromosome_name", "strand", 
# 		"external_gene_name", "cds_length", "genomic_coding_start", "genomic_coding_end", 
# 		"cds_start", "cds_end")
# cds_data = getBM(mart = mart, attributes = attr)
# cds_data[cds_data == ""] = NA # a few fields return empty (NA); need to drop these

# save data for personal use if frequently building hg19
# write.table(cds_data, "ensembl_hg19.txt", sep="\t", quote=F, row.names = F)

# pick 1,000 "random" genes 
# need to get all transcripts from each gene, so you can't just randomly sample rows
# set.seed(3019)
# all_genes = unique(cds_data$ensembl_gene_id)
# chosen = all_genes[sample(1:length(all_genes), 1000)]
# cds_small = cds_data[cds_data$ensembl_gene_id %in% chosen,]
# saveRDS(cds_small, "ensembl_cds_hg19_small.rds")


## Generate a mini RefCDS object for test purposes
cds_small = fread(get_test_file("cds_small_stop_codons_included.txt"))
refcds_and_gr_genes = build_RefCDS(cds_small, genome = "hg19", cds_ranges_lack_stop_codons = F)

saveRDS(refcds_and_gr_genes, "refcds_hg19_small.rds")
setwd(prev_dir)

