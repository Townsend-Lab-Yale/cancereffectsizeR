setwd(system.file("tests/testthat/testdata", package = "cancereffectsizeR"))


luad = load_maf(cesa = CESAnalysis(), maf = "luad.hg19.maf.txt", sample_col = "sample_id", tumor_allele_col = "Tumor_Seq_Allele2")
luad = calc_baseline_mutation_rates(luad, covariate_file = "lung_pca") 
saveRDS(luad, "cesa_for_snv.rds")

# save results to serve as expected test output
test_genes = c("TTN", "EGFR", "ASXL3", "KRAS", "RYR2", "USH2A", "CSMD3", "TP53", "CSMD1", "LRP1B", 
               "ZFHX4", "FAT3", "CNTNAP5", "PCDH15", "NEB", "RYR3", "DMD", "KATNAL1", 
               "OR13H1", "KSR1")
luad = effect_size_SNV(luad, genes = test_genes, analysis = "SNV")
results = selection_results_converter(luad)
saveRDS(results, "single_stage_snv_results.rds")


## Re-run, incorporating panel 1 test data
luad = load_maf(cesa = CESAnalysis(), maf = "luad.hg19.maf.txt", sample_col = "sample_id", tumor_allele_col = "Tumor_Seq_Allele2")

luad = load_maf(cesa = luad, maf = "panel1.hg19.maf.txt", 
	covered_regions = "tiny_panel.bed.txt")
luad = calc_baseline_mutation_rates(luad, covariate_file = "lung_pca") 

saveRDS(luad, "cesa_for_panel1_snv.rds")
test_genes = c("EGFR", "TP53", "KRAS", "OPTN", "ZFYVE27") # need all of theses genes for certain tests
luad =effect_size_SNV(luad, genes = test_genes)

results = selection_results_converter(luad)
saveRDS(results, "snv_with_panel_1.rds")

## Re-run again, incorporating panel 1 and panel 2 test data
luad = load_maf(cesa = CESAnalysis(), maf = "luad.hg19.maf.txt", sample_col = "sample_id", tumor_allele_col = "Tumor_Seq_Allele2")
luad = load_maf(cesa = luad, maf = "panel1.hg19.maf.txt", 
	covered_regions = "tiny_panel.bed.txt")
luad = load_maf(cesa = luad, maf = "panel2.hg19.maf.txt", 
	covered_regions = "tiny_panel.bed.txt")
luad = calc_baseline_mutation_rates(luad, covariate_file = "lung_pca")

saveRDS(luad, "cesa_for_panel_1_and_2_snv.rds")
luad = effect_size_SNV(luad, genes = test_genes)
results = selection_results_converter(luad)
saveRDS(results, "snv_with_panel_1_and_2.rds")



