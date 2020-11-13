prev_dir = setwd(system.file("tests/test_data/", package = "cancereffectsizeR"))

# read in the MAF, assign three random groups, and annotate
maf = fread("luad.hg19.maf.txt")
set.seed(879)
fruits = c("cherry", "marionberry", "mountain_apple")
maf[, group := sample(fruits, size = 1), by = "sample_id"]
luad = load_maf(cesa = CESAnalysis(ref_set = "ces_hg19_v1", sample_groups = fruits), maf = maf, group_col = "group", 
                sample_col = "sample_id")
saveRDS(luad, "annotated_fruit_cesa.rds")


# signatures and trinuc rates
luad = trinuc_mutation_rates(luad, cores = 4, signature_set = "COSMIC_v3.1",
                             signatures_to_remove = suggest_cosmic_signatures_to_remove("LUAD", treatment_naive = TRUE, quiet = TRUE))

fwrite(luad$trinuc_rates, "luad_hg19_trinuc_rates.txt", sep = "\t")
fwrite(luad$mutational_signatures$all, "luad_hg19_sig_table_with_artifacts.txt", sep = "\t")
fwrite(luad$mutational_signatures$biological, "luad_hg19_sig_table_biological.txt", sep = "\t")

# to generate test data for dndscv,
# run gene_mutation_rates(luad, covariates = "lung") with breakpoints before/after run_dndscv
# and save the input list and raw output to the .rds files

luad = gene_mutation_rates(luad, covariates = "lung", sample_group = "marionberry")
luad = gene_mutation_rates(luad, covariates = "lung", sample_group = c("cherry", "mountain_apple"))
fwrite(luad$gene_rates, "luad_fruit_gene_rates.txt", sep = "\t")


# save results to serve as expected test output
test_genes = c("EGFR", "ASXL3", "KRAS", "RYR2", "USH2A", "CSMD3", "TP53", "CSMD1", "LRP1B", 
               "ZFHX4", "FAT3", "CNTNAP5", "PCDH15", "NEB", "RYR3", "DMD", "KATNAL1", 
               "OR13H1", "KSR1")
luad = ces_variant(luad, variants = select_variants(luad, genes = test_genes))
fwrite(luad@selection_results[[1]], "fruit_sswm_out.txt", sep = "\t")

# Three big genes and a variant that is the only mutation in its gene in the data set
luad = ces_variant(luad, select_variants(luad, genes = c("EGFR", "KRAS", "TP53"), variant_passlist = "CR2 R247L"),
               model = "sswm_sequential", group_ordering = list(c("marionberry", "cherry"), "mountain_apple"))
fwrite(luad@selection_results[[2]], "fruit_sswm_sequential_out.txt", sep = "\t")

results = ces_gene_epistasis(luad, genes = c("EGFR", "KRAS", "TP53"), conf = .95)
saveRDS(results, "epistasis_results.rds")

setwd(prev_dir)






