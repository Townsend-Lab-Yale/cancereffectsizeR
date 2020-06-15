test_that("Trinucleotide signature weight calculation", {
  cesa = get_test_data("cesa_for_trinuc_weighting_calc.rds")
  to_remove = suggest_cosmic_v3_signatures_to_remove(cancer_type = "LUAD", treatment_naive = T, quiet = T)
  expect_identical(to_remove, c("SBS7a", "SBS7b", "SBS7c", "SBS7d", "SBS8", "SBS10a", "SBS10b", "SBS11", "SBS12", 
                                "SBS14", "SBS16", "SBS19", "SBS20", "SBS21", "SBS22", "SBS23", "SBS24", "SBS25", 
                                "SBS26", "SBS30", "SBS31", "SBS32", "SBS33", "SBS34", "SBS35", "SBS36", "SBS37", 
                                "SBS38", "SBS39", "SBS41", "SBS42", "SBS44"))
  cesa = trinuc_mutation_rates(cesa, signatures_to_remove = to_remove)
  trinuc_ak = get_test_data("trinuc_mut_weighting.rds")
  expect_equal(cesa@trinucleotide_mutation_weights, trinuc_ak)
})

test_that("dNdScv and MAF annotation", {
  cesa = get_test_data("cesa_for_dndscv_and_anno.rds")
  dndscv_input = cancereffectsizeR:::dndscv_preprocess(cesa = cesa, covariate_file = "lung_pca")
  dndscv_input_ak = get_test_data("dndscv_input_single.rds")
  expect_equal(dndscv_input[[1]][-4], dndscv_input_ak[[1]][-4]) # refdb path will very on dev/prod due to inst dir
  dndscv_output = get_test_data("dndscv_raw_output_single.rds")
  cesa = cancereffectsizeR:::dndscv_postprocess(cesa, dndscv_output)
  sel_cv = get_test_data("sel_cv.rds")
  expect_equal(cesa@dndscv_out_list[[1]]$sel_cv, sel_cv)
  mutrates = get_test_data("mutrates.rds")
  expect_equal(cesa@mutrates, mutrates)
  cesa = annotate_variants(cesa)
  annotated_maf = get_test_data("annotated_maf_df.rds")
  expect_equal(cesa@maf, annotated_maf)
})

# will use SNV-analysis-ready object for remaining tests
cesa = get_test_data("cesa_for_snv.rds")

test_that("Handle missing or invalid gene choice in SNV analysis", {
  # Error when any requested gene is not in RefCDS data
  expect_error(ces_snv(cesa, genes = c("KRAS", "TP53", "notagene")),
               "requested genes were not found")

  # AC006486.1 is not in the data set; error when no genes requested are in the data set
  expect_error(ces_snv(cesa, genes = c("AC006486.1")),
               "None of the requested genes have eligible mutations")

  # Expect a message when one or more of the genes requested isn't in data set
  # This call quits early after receiving the message to save time
  expect_match(tryCatch({ces_snv(cesa, genes = c("AC006486.1", "TP53"))},
                        message = function(m) {m$message}),
               "The following requested genes have no eligible mutations")
})


# genes to plug into ces_snv; some with high-effect-size SNVs, others random
test_genes = c("TTN", "EGFR", "ASXL3", "KRAS", "RYR2", "USH2A", "CSMD3", "TP53", "CSMD1", "LRP1B", 
               "ZFHX4", "FAT3", "CNTNAP5", "PCDH15", "NEB", "RYR3", "DMD", "KATNAL1", 
               "OR13H1", "KSR1")
test_that("SNV effect size calculation", {
  cesa = ces_snv(cesa, genes = test_genes, include_nonrecurrent_variants = T)
  results = cesa@selection_results
  results_ak = get_test_data("single_stage_snv_results.rds")
  
  # selection results will vary slightly by machine architecture (probably due to numerical precision issues)
  expect_equal(results, results_ak, tolerance = 1e-5)
})


test_that("Gene-level SNV epistasis analysis", {
  cesa = ces_gene_epistasis(cesa, genes = c("EGFR", "KRAS", "TP53"), return_all_opm_output = T)
  results = cesa@gene_epistasis_results
  results_ak = get_test_data("epistasis_results.rds")
  expect_equal(results, results_ak, tolerance = 1e-4)
  opm_ak = get_test_data("epistasis_opm.rds")
  expect_equal(cesa@advanced$opm_output[,!"xtime"], opm_ak[,!"xtime"], tolerance = 1e-4) # runtime will vary
})






