test_that("Trinucleotide signature weight calculation", {
  cesa = get_test_data("cesa_for_trinuc_weighting_calc.rds")
  cesa = trinucleotide_mutation_weights(cesa)
  trinuc_ak = get_test_data("trinuc_mut_weighting.rds")
  expect_equal(cesa@trinucleotide_mutation_weights, trinuc_ak)
})

test_that("dNdScv and MAF annotation", {
  cesa = get_test_data("cesa_for_dndscv_and_anno.rds")
  cesa = expect_warning(gene_level_mutation_rates(cesa, covariate_file = "lung_pca"), "Same mutations observed in different sample")
  sel_cv = get_test_data("sel_cv.rds")
  expect_equal(data.table(cesa@dndscv_out_list$`1`$sel_cv), data.table(sel_cv))
  mutrates = get_test_data("mutrates.rds")
  expect_equal(cesa@mutrates_list$`1`, mutrates)
  cesa = annotate_gene_maf(cesa)
  annotated_maf = get_test_data("annotated_maf_df.rds")
  expect_equal(data.frame(cesa@annotated.snv.maf), data.frame(annotated_maf))
})

# will use SNV-analysis-ready object for remaining tests
cesa = get_test_data("cesa_for_snv.rds")

test_that("Handle missing or invalid gene choice in SNV analysis", {
  # Error when any requested gene is not in RefCDS data
  expect_error(ces_snv(cesa, genes = c("KRAS", "TP53", "notagene")),
               "requested genes were not found")
  
  # AC006486.1 is not in the data set; error when no genes requested are in the data set
  expect_error(ces_snv(cesa, genes = c("AC006486.1")),
               "None of the requested genes have mutations")
  
  # Expect a message when one or more of the genes requested isn't in data set
  # This call quits early after receiving the message to save time
  expect_match(tryCatch({ces_snv(cesa, genes = c("AC006486.1", "TP53"))}, 
                        message = function(m) {m$message}),
               "The following requested genes have no mutations")
})


# genes to plug into ces_snv; some with high-effect-size SNVs, others random
test_genes = c("TTN", "EGFR", "ASXL3", "KRAS", "RYR2", "USH2A", "CSMD3", "TP53", "CSMD1", "LRP1B", 
               "ZFHX4", "FAT3", "CNTNAP5", "PCDH15", "NEB", "RYR3", "DMD", "KATNAL1", 
               "OR13H1", "KSR1")
test_that("SNV effect size calculation", {
  cesa = ces_snv(cesa, genes = test_genes)
  results = selection_results_converter(cesa)
  results_ak = get_test_data("single_stage_snv_results.rds")
  
  # selection results will vary slightly by machine due to numerical precision issues during parameter optimization
  ## For now, need to sort dfs by column name because dplyr version affects column orderings
  expect_equal(results[,order(names(results))], results_ak[,order(names(results_ak))], tolerance = 1e-5)
})




test_that("Gene-level SNV epistasis analysis", {
  cesa = ces_gene_epistasis(cesa, genes = c("EGFR", "KRAS", "TP53"))
  results = cesa@selection_results
  results_ak = get_test_data("epistasis_results.rds")
  expect_equal(results, results_ak, tolerance = 1e-7)
})






