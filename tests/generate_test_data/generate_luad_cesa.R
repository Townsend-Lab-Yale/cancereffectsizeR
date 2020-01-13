setwd(system.file("tests/testthat/testdata", package = "cancereffectsizeR"))

# read in the MAF used for all the testing
luad = load_maf(cesa = CESAnalysis(genome="hg19"), maf = "luad.hg19.maf.txt", sample_col = "sample_id", tumor_allele_col = "Tumor_Seq_Allele2")

# run all pre-processing and save for quick import
luad = calc_baseline_mutation_rates(luad, covariate_file = "lung_pca")
saveRDS(luad, "cesa_for_snv.rds")

# save results to serve as expected test output
test_genes = c("TTN", "EGFR", "ASXL3", "KRAS", "RYR2", "USH2A", "CSMD3", "TP53", "CSMD1", "LRP1B", 
               "ZFHX4", "FAT3", "CNTNAP5", "PCDH15", "NEB", "RYR3", "DMD", "KATNAL1", 
               "OR13H1", "KSR1")
luad = ces_snv(luad, genes = test_genes)
results = selection_results_converter(luad)
saveRDS(results, "single_stage_snv_results.rds")
luad = ces_gene_epistasis(luad, genes = c("EGFR", "KRAS", "TP53"))
saveRDS(luad@selection_results, "epistasis_results.rds")


# repeat with smaller data set for dndscv testing
dndscv_samples = c("sample-1", "sample-106", "sample-108", "sample-11", "sample-31", "sample-33", "sample-35", 
                   "sample-40", "sample-46", "sample-6", "sample-67", "sample-68", "sample-7", "sample-71", 
                   "sample-73", "sample-74", "sample-77", "sample-82", "sample-83", "sample-95", "sample-99")

maf_for_dndscv = luad@maf[Unique_Patient_Identifier %in% dndscv_samples]

for_dndscv = load_maf(cesa = CESAnalysis(genome="hg19"), maf = maf_for_dndscv, sample_col = "Unique_Patient_Identifier", tumor_allele_col = "Tumor_Allele")
for_dndscv = trinucleotide_mutation_weights(for_dndscv)
saveRDS(for_dndscv, "cesa_for_dndscv_and_anno.rds")

dndscv_out = gene_level_mutation_rates(for_dndscv, covariate_file = "lung_pca")
saveRDS(dndscv_out@dndscv_out_list$`1`$sel_cv, "sel_cv.rds")
saveRDS(dndscv_out@mutrates_list$`1`, "mutrates.rds")
anno_out = annotate_gene_maf(dndscv_out)
saveRDS(anno_out@annotated.snv.maf, "annotated_maf_df.rds")


# ## Re-run, incorporating panel 1 test data
# luad = load_maf(cesa = CESAnalysis(), maf = "luad.hg19.maf.txt", sample_col = "sample_id", tumor_allele_col = "Tumor_Seq_Allele2")
# 
# luad = load_maf(cesa = luad, maf = "panel1.hg19.maf.txt", 
# 	covered_regions = "tiny_panel.bed.txt")
# luad = calc_baseline_mutation_rates(luad, covariate_file = "lung_pca") 
# 
# saveRDS(luad, "cesa_for_panel1_snv.rds")
# test_genes = c("EGFR", "TP53", "KRAS", "OPTN", "ZFYVE27") # need all of theses genes for certain tests
# luad =ces_snv(luad, genes = test_genes)
# 
# results = selection_results_converter(luad)
# saveRDS(results, "snv_with_panel_1.rds")
# 
# ## Re-run again, incorporating panel 1 and panel 2 test data
# luad = load_maf(cesa = CESAnalysis(), maf = "luad.hg19.maf.txt", sample_col = "sample_id", tumor_allele_col = "Tumor_Seq_Allele2")
# luad = load_maf(cesa = luad, maf = "panel1.hg19.maf.txt", 
# 	covered_regions = "tiny_panel.bed.txt")
# luad = load_maf(cesa = luad, maf = "panel2.hg19.maf.txt", 
# 	covered_regions = "tiny_panel.bed.txt")
# luad = calc_baseline_mutation_rates(luad, covariate_file = "lung_pca")
# 
# saveRDS(luad, "cesa_for_panel_1_and_2_snv.rds")
# luad = ces_snv(luad, genes = test_genes)
# results = selection_results_converter(luad)
# saveRDS(results, "snv_with_panel_1_and_2.rds")



