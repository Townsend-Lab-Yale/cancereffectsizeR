cesa = load_cesa(get_test_file("cesa_for_trinuc_weighting_calc.rds"))
test_that("Old signature set validates", 
          expect_no_error(validate_signature_set(ces.refset.hg19$signatures$COSMIC_v3.1)))
test_that("Trinucleotide signature weight calculation", {
  to_remove = suggest_cosmic_signature_exclusions(cancer_type = "LUAD", treatment_naive = T, quiet = T)
  expect_identical(to_remove, c("SBS7a", "SBS7b", "SBS7c", "SBS7d", "SBS8", "SBS10a", "SBS10b", "SBS11", "SBS12",
                                "SBS14", "SBS16", "SBS19", "SBS20", "SBS21", "SBS22", "SBS22a", "SBS22b", "SBS23", "SBS24", "SBS25",
                                "SBS26", "SBS30", "SBS31", "SBS32", "SBS33", "SBS34", "SBS35", "SBS36", "SBS37",
                                "SBS38", "SBS39", "SBS41", "SBS42", "SBS44", "SBS84", "SBS85", "SBS86", "SBS87",
                                "SBS88", "SBS89", "SBS90", "SBS91", "SBS92", "SBS93", "SBS94", "SBS99"))
  
  cesa = trinuc_mutation_rates(cesa, signature_exclusions = to_remove, signature_set = "COSMIC_v3.1",
                               signature_extractor = 'MutationalPatterns')
  
  expect_error(trinuc_mutation_rates(cesa, signature_set = "COSMIC_v3.1"), 'have already been run')
  expect_error(trinuc_mutation_rates(cesa, signature_set = "COSMIC_v3.1", samples = cesa$samples[1]), 'have already been run')
  
  trinuc_ak = get_test_data("trinuc_mut_weighting.rds")
  expect_equal(length(waldo::compare(cesa@trinucleotide_mutation_weights, trinuc_ak, 
                                     tolerance = 1e-4, ignore_attr = 'index')), 0)
  
  
  # Ensure SNV counts (total and used by MP) look right
  expect_identical(cesa@trinucleotide_mutation_weights$signature_weight_table[, unlist(.(total_snvs, sig_extraction_snvs))],
                   c(67, 181, 38, 18, 0, 66, 180, 38, 18, 0))
  
  
  full_weight_table = get_signature_weights(cesa)
  expect_equal(full_weight_table[, .N], 5)
  expect_equal(full_weight_table[Unique_Patient_Identifier == "one_indel", unlist(.(total_snvs, sig_extraction_snvs))], c(0, 0))

  # All tumors should have an above-threshold number of usable SNVs, or they should be group-average-blended (but never both)
  expect_true(all(full_weight_table[, xor(sig_extraction_snvs > 49, group_avg_blended == T)]))
  
  # test set_trinuc_rates and set_signature_weights
  trinuc_rates = cesa$trinuc_rates
  prev_rates_matrix = cesa@trinucleotide_mutation_weights$trinuc_proportion_matrix
  prev_relative_bio_sig = cesa@trinucleotide_mutation_weights$signature_weight_table
  prev_raw_sig = cesa@trinucleotide_mutation_weights$signature_weight_table_with_artifacts
  cesa = clear_trinuc_rates_and_signatures(cesa)
  
  
  
  cesa2 = set_trinuc_rates(cesa, trinuc_rates = trinuc_rates)
  expect_equal(prev_rates_matrix, cesa2@trinucleotide_mutation_weights$trinuc_proportion_matrix)
  
  # setting signature weights using raw weights will reproduce the same 
  # relative biological weights and trinuc_rates for non-mean-blended samples
  raw_weight_input = prev_raw_sig[, .SD, .SDcols = patterns("(SBS)|(Uni)")]
  summed_weight = rowSums(raw_weight_input[, -"Unique_Patient_Identifier"])
  raw_weight_input[, names(.SD) := lapply(.SD, '/', summed_weight), .SDcols = patterns('SBS')]
  cesa3 = set_signature_weights(cesa, signature_set = "COSMIC_v3.1", weights = raw_weight_input)
  
  
  which_pure = prev_relative_bio_sig[group_avg_blended == F, which = T]
  expect_equal(cesa3@trinucleotide_mutation_weights$signature_weight_table[which_pure, -c("sig_extraction_snvs", "group_avg_blended")],
            prev_relative_bio_sig[which_pure, -c("sig_extraction_snvs", "group_avg_blended")])
  expect_equivalent(cesa3@trinucleotide_mutation_weights$raw_signature_weights[,  .SD, .SDcols = patterns("(SBS)|(Uni)")],
            raw_weight_input)
  expect_equal(cesa3@trinucleotide_mutation_weights$trinuc_proportion_matrix[which_pure,],
            prev_rates_matrix[which_pure, ])
  
  # setting with the previous biological weights will yield identical weights/rates to the previous run
  cesa4 = set_signature_weights(cesa, signature_set = "COSMIC_v3.1", weights = prev_relative_bio_sig)
  expect_equal(cesa4@trinucleotide_mutation_weights$signature_weight_table[, -c("sig_extraction_snvs", "group_avg_blended")],
            prev_relative_bio_sig[, -c("sig_extraction_snvs", "group_avg_blended")])
  
  current = cesa4@trinucleotide_mutation_weights$trinuc_proportion_matrix
  expect_equal(current[sort(rownames(current)),], prev_rates_matrix[sort(rownames(prev_rates_matrix)),], 
               tolerance = 1e-4) # weird that it's not a little more precise, but whatever
  
  
  # All tumors should have an above-threshold number of usable SNVs, or they should be group-average-blended (but never both)
  expect_true(all(full_weight_table[, xor(sig_extraction_snvs > 49, group_avg_blended == T)]))
})

test_that('trinuc_snv_counts output is usable by extractors', {
  mp_counts = trinuc_snv_counts(cesa$maf, genome = cesa$reference_data$genome, style = 'MutationalPatterns')
  expect_silent(MutationalPatterns::fit_to_signatures(mut_matrix = mp_counts, signatures = t(ces.refset.hg19$signatures$COSMIC_v3$signatures)))
})


