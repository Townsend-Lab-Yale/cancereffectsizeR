luad = CESAnalysis(progression_order = 1:4)
luad = load_maf(luad, maf = "tests/testthat/testdata/luad.hg19.maf.txt", sample_col = "sample_id", 
                tumor_allele_col = "Tumor_Seq_Allele2", progression_col = "fake_stage")
luad = calc_baseline_mutation_rates(luad, covariate_file = "lung_pca")
saveRDS(luad, "tests/testthat/testdata/cesa_for_snv_multi.rds")
