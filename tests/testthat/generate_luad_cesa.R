luad = CESAnalysis(maf = "tests/testthat/testdata/luad.hg19.maf.txt", sample_col = "sample_id", tumor_allele_col = "Tumor_Seq_Allele2")
luad = calc_baseline_mutation_rates(luad, covariate_file = "lung_pca", BIG_mem = F) # BIG_mem  option will likely be removed in future with its effects made default
saveRDS(luad, "cesa_for_snv.rds")
