prev_dir = setwd(system.file("tests/test_data/", package = "cancereffectsizeR"))

# read in the MAF, assign three random groups, and annotate

maf = fread("brca_subset_hg38.maf.gz")
maf = preload_maf(maf, refset = 'ces.refset.hg38')

sample_data = data.table(Unique_Patient_Identifier = unique(maf$Unique_Patient_Identifier))
set.seed(879)
fruits = c("cherry", "marionberry", "mountain_apple")
sample_data[, fruit := sample(fruits, .N, replace = T)]

cesa = load_maf(cesa = CESAnalysis(refset = "ces.refset.hg38"), maf = maf, maf_name = "BRCA")
cesa = load_sample_data(cesa, sample_data)

tgs_maf_file = "../../inst/tutorial/metastatic_breast_2021_hg38.maf"
tgs_maf = preload_maf(maf = tgs_maf_file, refset = "ces.refset.hg38")
tgs_maf[, fruit := sample(fruits, 1), by = 'Unique_Patient_Identifier']


# Define coverage using the coding regions of these genes, as defined by the refset
top_tgs_genes <- c(
  "TP53", "PIK3CA", "ESR1", "CDH1", "GATA3", "KMT2C", 'KRAS',
  "MAP3K1", "AKT1", "ARID1A", "FOXA1", "TBX3", "PTEN", 'EGFR'
)
tgs_coverage <- ces.refset.hg38$gr_genes[ces.refset.hg38$gr_genes$gene %in% top_tgs_genes]
tgs_maf$pM <- "M1"
cesa <- load_maf(cesa,
                 maf = tgs_maf, sample_data_cols = c("pM", 'fruit'), maf_name = "MBC",
                 coverage = "targeted", covered_regions = tgs_coverage,
                 covered_regions_name = "top_genes", covered_regions_padding = 10
)

# signatures and trinuc rates
cesa = trinuc_mutation_rates(cesa = cesa, signature_set = "COSMIC_v3.2",
                             signature_exclusions = suggest_cosmic_signature_exclusions("LUAD", treatment_naive = TRUE, quiet = TRUE))

# fwrite(luad$trinuc_rates, "luad_hg19_trinuc_rates.txt", sep = "\t")
# fwrite(luad$mutational_signatures$raw_attributions, "luad_hg19_sig_table_with_artifacts.txt", sep = "\t")
# fwrite(luad$mutational_signatures$biological_weights, "luad_hg19_sig_table_biological.txt", sep = "\t")

cesa = gene_mutation_rates(cesa, covariates = "lung", samples = cesa$samples[fruit == "marionberry"])
cesa = gene_mutation_rates(cesa, covariates = "lung", samples = cesa$samples[fruit %in% c("cherry", "mountain_apple")])

save_cesa(cesa, 'cesa_hg38_for_test.rds')

# save results to serve as expected test output
test_genes = c("EGFR", "ASXL3", "KRAS", "RYR2", "USH2A", "CSMD3", "TP53", 
               "CSMD1", "LRP1B", "ZFHX4", "FAT3", "CNTNAP5", "PCDH15", "NEB", 
               "RYR3", "PIK3CA", "ESR1", 
               "MAP3K1", "AKT1", "ARID1A", "FOXA1", 
               "TBX3", "PTEN")
cesa = ces_variant(cesa, variants = select_variants(cesa, genes = test_genes))
fwrite(cesa@selection_results[[1]], "default_model_effects_brca_hg38.txt", sep = "\t")




variants_to_use = cesa$variants[gene %in% c('EGFR', 'KRAS', 'TP53') & samples_covering == cesa$samples[, .N]]
cesa = ces_gene_epistasis(cesa, genes = c("EGFR", "KRAS", "TP53"), variants = variants_to_use, conf = .95)
saveRDS(cesa$epistasis[[1]], "epistatic_effects.rds")

setwd(prev_dir)






