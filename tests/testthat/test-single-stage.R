cesa = load_cesa(get_test_file("cesa_for_trinuc_weighting_calc.rds"))
test_that("Trinucleotide signature weight calculation", {
  to_remove = suggest_cosmic_signatures_to_remove(cancer_type = "LUAD", treatment_naive = T, quiet = T)
  expect_identical(to_remove, c("SBS7a", "SBS7b", "SBS7c", "SBS7d", "SBS8", "SBS10a", "SBS10b", "SBS11", "SBS12", 
                                "SBS14", "SBS16", "SBS19", "SBS20", "SBS21", "SBS22", "SBS23", "SBS24", "SBS25", 
                                "SBS26", "SBS30", "SBS31", "SBS32", "SBS33", "SBS34", "SBS35", "SBS36", "SBS37", 
                                "SBS38", "SBS39", "SBS41", "SBS42", "SBS44", "SBS84", "SBS85", "SBS86", "SBS87", 
                                "SBS88", "SBS90"))
  
  cesa = trinuc_mutation_rates(cesa, signatures_to_remove = to_remove, signature_set = "COSMIC_v3.1")
  
  # Ensure SNV counts (total and used by dS) look right
  expect_identical(cesa@trinucleotide_mutation_weights$signature_weight_table[, unlist(.(total_snvs, sig_extraction_snvs))],
                   c(67, 181, 38, 66, 180, 38))
  
  
  full_weight_table = get_signature_weights(cesa, include_tumors_without_data = T)
  expect_equal(full_weight_table[, .N], 5)
  expect_equal(full_weight_table[Unique_Patient_Identifier == "one_indel", unlist(.(total_snvs, sig_extraction_snvs))], c(0, 0))
  expect_equal(full_weight_table[Unique_Patient_Identifier == "zeroed-out", unlist(.(total_snvs, sig_extraction_snvs))], c(18, 0))
  
  # All tumors should have an above-threshold number of usable SNVs, or they should be group-average-blended (but never both)
  expect_true(all(full_weight_table[, xor(sig_extraction_snvs > 49, group_avg_blended == T)]))
  
  
  trinuc_ak = get_test_data("trinuc_mut_weighting.rds")
  expect_equal(cesa@trinucleotide_mutation_weights, trinuc_ak)
  
  trinuc_rates = cesa$trinuc_rates
  prev_rates_matrix = cesa@trinucleotide_mutation_weights$trinuc_proportion_matrix
  cesa@trinucleotide_mutation_weights = list()
  
  cesa2 = set_trinuc_rates(cesa, trinuc_rates = trinuc_rates)
  expect_equal(prev_rates_matrix, cesa2@trinucleotide_mutation_weights$trinuc_proportion_matrix)
  
  
})

# will use SNV-analysis-ready object for remaining tests
cesa = load_cesa(get_test_file("cesa_for_snv.rds"))

test_that("gene mutation rates", {
  dndscv_input = cancereffectsizeR:::dndscv_preprocess(cesa = cesa, covariates = "lung")
  dndscv_input_ak = get_test_data("dndscv_input_single.rds")
  
  # all but first five columns of mutations ignored by dndscv (later, should get rid of those columns in gene_mutation_rates)
  dndscv_input[[1]]$mutations = dndscv_input[[1]]$mutations[, 1:5]
  dndscv_input_ak[[1]]$mutations = dndscv_input_ak[[1]]$mutations[, 1:5]
  expect_equal(dndscv_input[[1]][-4], dndscv_input_ak[[1]][-4]) # refdb path will very on dev/prod due to inst dir
  dndscv_output = get_test_data("dndscv_raw_output_single.rds")
  dndscv_output = cancereffectsizeR:::dndscv_postprocess(cesa, dndscv_output)
  sel_cv = get_test_data("sel_cv.rds")
  expect_equal(dndscv_output@dndscv_out_list[[1]], sel_cv)
  mutrates = get_test_data("mutrates.rds")
  expect_equal(dndscv_output@mutrates, mutrates)
})


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
test_genes = c("EGFR", "ASXL3", "KRAS", "RYR2", "USH2A", "CSMD3", "TP53", "CSMD1", "LRP1B", 
               "ZFHX4", "FAT3", "CNTNAP5", "PCDH15", "NEB", "RYR3", "DMD", "KATNAL1", 
               "OR13H1", "KSR1")
test_that("SNV effect size calculation", {
  cesa = ces_snv(cesa, genes = test_genes, include_nonrecurrent_variants = T, cores = 1)
  results = cesa@selection_results
  results_ak = get_test_data("single_stage_snv_results.rds")
  
  # selection results will vary slightly by machine architecture (probably due to numerical precision issues)
  expect_equal(results, results_ak, tolerance = 1e-5, ignore.row.order = T)
})

test_that("ces_snv with user-supplied variants", {
  expect_error(ces_snv(cesa, variant = list(snv_id = "10:100190376_C>A")), "No variants pass filters")
  expect_equal(ces_snv(cesa, variant = list(snv_id = "10:100190376_C>A"), include_nonrecurrent_variants = T)@selection_results[, .N], 1)
  expect_error(ces_snv(cesa, genes = "TP53", variant = list(aac_id = "KRAS_G12D_ENSP00000256078"), "No variants pass filters"))
  expect_equal(ces_snv(cesa, genes = c("TP53", "KRAS"), variant = list(aac_id = "KRAS_G12D_ENSP00000256078"))@selection_results[, .N], 1)
})

test_that("Gene-level SNV epistasis analysis", {
  cesa = ces_gene_epistasis(cesa, genes = c("EGFR", "KRAS", "TP53"))
  results = cesa@gene_epistasis_results
  results_ak = get_test_data("epistasis_results.rds")
  expect_equal(results, results_ak, tolerance = 1e-3)
})






