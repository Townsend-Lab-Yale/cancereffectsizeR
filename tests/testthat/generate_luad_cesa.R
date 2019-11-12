luad = CESAnalysis(maf = "tests/testthat/testdata/luad.hg19.maf.txt", sample_col = "sample_id", tumor_allele_col = "Tumor_Seq_Allele2")
luad = calc_baseline_mutation_rates(luad, covariate_file = "lung_pca") 
saveRDS(luad, "cesa_for_snv.rds")
