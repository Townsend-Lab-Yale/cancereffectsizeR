cesa = load_cesa(get_test_file("cesa_for_trinuc_weighting_calc.rds"))
test_that("Trinucleotide signature weight calculation", {
  to_remove = suggest_cosmic_signatures_to_remove(cancer_type = "LUAD", treatment_naive = T, quiet = T)
  expect_identical(to_remove, c("SBS7a", "SBS7b", "SBS7c", "SBS7d", "SBS8", "SBS10a", "SBS10b", "SBS11", "SBS12",
                                "SBS14", "SBS16", "SBS19", "SBS20", "SBS21", "SBS22", "SBS23", "SBS24", "SBS25",
                                "SBS26", "SBS30", "SBS31", "SBS32", "SBS33", "SBS34", "SBS35", "SBS36", "SBS37",
                                "SBS38", "SBS39", "SBS41", "SBS42", "SBS44", "SBS84", "SBS85", "SBS86", "SBS87",
                                "SBS88", "SBS89", "SBS90"))
  
  cesa = trinuc_mutation_rates(cesa, signatures_to_remove = to_remove, signature_set = "COSMIC_v3.1")
  trinuc_ak = get_test_data("trinuc_mut_weighting.rds")
  expect_equal(cesa@trinucleotide_mutation_weights, trinuc_ak)
  
  
  
  # Ensure SNV counts (total and used by dS) look right
  expect_identical(cesa@trinucleotide_mutation_weights$signature_weight_table[, unlist(.(total_snvs, sig_extraction_snvs))],
                   c(67, 181, 38, 18, 0, 66, 180, 38, 0, 0))
  
  
  full_weight_table = get_signature_weights(cesa)
  expect_equal(full_weight_table[, .N], 5)
  expect_equal(full_weight_table[Unique_Patient_Identifier == "one_indel", unlist(.(total_snvs, sig_extraction_snvs))], c(0, 0))
  expect_equal(full_weight_table[Unique_Patient_Identifier == "zeroed-out", unlist(.(total_snvs, sig_extraction_snvs))], c(18, 0))
  
  # All tumors should have an above-threshold number of usable SNVs, or they should be group-average-blended (but never both)
  expect_true(all(full_weight_table[, xor(sig_extraction_snvs > 49, group_avg_blended == T)]))
  
  # test set_trinuc_rates and set_signature_weights
  trinuc_rates = cesa$trinuc_rates
  prev_rates_matrix = cesa@trinucleotide_mutation_weights$trinuc_proportion_matrix
  prev_relative_bio_sig = cesa@trinucleotide_mutation_weights$signature_weight_table
  prev_raw_sig = cesa@trinucleotide_mutation_weights$signature_weight_table_with_artifacts
  cesa@trinucleotide_mutation_weights = list()
  
  cesa2 = set_trinuc_rates(cesa, trinuc_rates = trinuc_rates)
  expect_equal(prev_rates_matrix, cesa2@trinucleotide_mutation_weights$trinuc_proportion_matrix)
  
  # setting signature weights using raw weights will reproduce the same 
  # relative biological weights and trinuc_rates for non-mean-blended samples
  cesa@trinucleotide_mutation_weights = list()
  raw_weight_input = prev_raw_sig[, .SD, .SDcols = patterns("(SBS)|(Uni)")]
  cesa3 = set_signature_weights(cesa, signature_set = "COSMIC_v3.1", weights = raw_weight_input)
  
  
  which_pure = prev_relative_bio_sig[group_avg_blended == F, which = T]
  all.equal(cesa3@trinucleotide_mutation_weights$signature_weight_table[which_pure, -c("sig_extraction_snvs", "group_avg_blended")],
            prev_relative_bio_sig[which_pure, -c("sig_extraction_snvs", "group_avg_blended")])
  all.equal(cesa3@trinucleotide_mutation_weights$signature_weight_table_with_artifacts[,  .SD, .SDcols = patterns("(SBS)|(Uni)")],
            raw_weight_input)
  all.equal(cesa3@trinucleotide_mutation_weights$trinuc_proportion_matrix[which_pure,],
            prev_rates_matrix[which_pure, ])
  
  # setting with the biological weights should get the same signature weights results,
  # but due to quirks in trinuc rate method, the trinuc rates of blended tumors will change slightly
  cesa4 = set_signature_weights(cesa, signature_set = "COSMIC_v3.1", weights = prev_relative_bio_sig)
  expect_equal(cesa4@trinucleotide_mutation_weights$signature_weight_table[, -c("sig_extraction_snvs", "group_avg_blended")],
            prev_relative_bio_sig[, -c("sig_extraction_snvs", "group_avg_blended")])
  
  expect_equal(cesa4@trinucleotide_mutation_weights$trinuc_proportion_matrix[which_pure, ], 
               prev_rates_matrix[which_pure,])
})