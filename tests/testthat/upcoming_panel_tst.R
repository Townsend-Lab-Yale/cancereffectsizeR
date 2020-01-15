test_genes = c("EGFR", "TP53", "KRAS", "OPTN", "ZFYVE27") # need all of theses genes for certain tests
luad = ces_snv(luad, genes = test_genes)

results = selection_results_converter(luad)
saveRDS(results, "snv_with_panel_1.rds")
q
## Re-run again, incorporating panel 1 and panel 2 test data
luad = load_maf(cesa = CESAnalysis(), maf = "luad.hg19.maf.txt", sample_col = "sample_id", tumor_allele_col = "Tumor_Seq_Allele2")
luad = load_maf(cesa = luad, maf = "panel1.hg19.maf.txt", 
                covered_regions = "tiny_panel.bed.txt")
luad = load_maf(cesa = luad, maf = "panel2.hg19.maf.txt", 
                covered_regions = "tiny_panel.bed.txt")
luad = calc_baseline_mutation_rates(luad, covariate_file = "lung_pca")

saveRDS(luad, "cesa_for_panel_1_and_2_snv.rds")
luad = ces_snv(luad, genes = test_genes)
results = selection_results_converter(luad)
saveRDS(results, "snv_with_panel_1_and_2.rds")