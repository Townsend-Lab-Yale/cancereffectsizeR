prev_dir = setwd(system.file("tests/test_data/", package = "cancereffectsizeR"))


# Other testing uses hg38; this still uses hg19 to check that ces.refset.hg19 is working
maf = preload_maf(maf = "luad.hg19.maf.txt", refset = 'ces.refset.hg19', sample_col = 'sample_id')
luad = load_maf(cesa = CESAnalysis(refset = "ces.refset.hg19"), maf = maf)
saveRDS(luad, "luad_hg19_for_gene_rate_test.rds")

# To generate test data for test-gene-rates.R,
# run the gene_mutation_rates call below with breakpoints before/after run_dndscv
# and save the input list and raw output to the .rds files. Also save the fit object.

luad = gene_mutation_rates(luad, covariates = "lung")
# saveRDS(dndscv_args[c('mutations', 'gene_list', 'cv')], 'dndscv_input_single.rds') # subset of run_dndscv() input.
# saveRDS(dndscv_output, "dndscv_raw_output_single.rds") # output of run_dndscv().
# saveRDS(fit, 'dndscv_test_fit_object.rds') # from get_dndscv_model_fit() in gene_mutation_rates()

saveRDS(luad@mutrates, "gene_mutation_rates_with_ci.rds")
setwd(prev_dir)






