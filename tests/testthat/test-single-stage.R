source(system.file("tests/testthat/helpers.R", package = "cancereffectsizeR"))
test_that("Trinucleotide signature weight calculation", {
  cesa = get_test_data("cesa_for_trinuc_weighting_calc.rds")
  cesa = expect_warning(trinucleotide_mutation_weights(cesa), "Some samples have fewer than 50 mutations")
  trinuc_ak = get_test_data("trinuc_mut_weighting.rds")
  expect_identical(cesa@trinucleotide_mutation_weights, trinuc_ak)
})

test_that("dNdScv and MAF annotation", {
  cesa = get_test_data("cesa_for_dndscv_and_anno.rds")
  cesa = expect_warning(gene_level_mutation_rates(cesa, covariate_file = "lung_pca"), "Same mutations observed in different sample")
  sel_cv = get_test_data("sel_cv.rds")
  expect_identical(cesa@dndscv_out_list$`1`$sel_cv, sel_cv)
  mutrates = get_test_data("mutrates.rds")
  expect_identical(cesa@mutrates_list$`1`, mutrates)
  cesa = annotate_gene_maf(cesa)
  annotated_maf = get_test_data("annotated_maf_df.rds")
  expect_identical(data.frame(cesa@annotated.snv.maf), annotated_maf)
})


# genes to plug into effect_size_SNV; some likely to be high effect size, others random
test_genes = c("TTN", "KRAS", "RYR2", "USH2A", "CSMD3", "TP53", "CSMD1", "LRP1B", 
               "ZFHX4", "FAT3", "CNTNAP5", "PCDH15", "NEB", "RYR3", "DMD", "KATNAL1", 
               "OR13H1", "KSR1", "VWA5B1", "IFNAR1", "ARRB2", "CMSS1", "SLC10A7", "ENOX2", "IFITM2")
test_that("SNV effect size calculation", {
  cesa = get_test_data("cesa_for_snv.rds")
  cesa = effect_size_SNV(cesa, genes = test_genes, analysis = "SNV", cores = 1)
  results = selection_results_converter(cesa)
  results_ak = get_test_data("single_stage_snv_results.rds")
  expect_identical(results, results_ak)
})




