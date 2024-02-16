cesa = load_cesa(get_test_file("cesa_for_trinuc_weighting_calc.rds"))
test_that("Trinucleotide signature weight calculation", {
  to_remove = suggest_cosmic_signature_exclusions(cancer_type = "LUAD", treatment_naive = T, quiet = T)
  expect_identical(to_remove, c("SBS7a", "SBS7b", "SBS7c", "SBS7d", "SBS8", "SBS10a", "SBS10b", "SBS11", "SBS12",
                                "SBS14", "SBS16", "SBS19", "SBS20", "SBS21", "SBS22", "SBS23", "SBS24", "SBS25",
                                "SBS26", "SBS30", "SBS31", "SBS32", "SBS33", "SBS34", "SBS35", "SBS36", "SBS37",
                                "SBS38", "SBS39", "SBS41", "SBS42", "SBS44", "SBS84", "SBS85", "SBS86", "SBS87",
                                "SBS88", "SBS89", "SBS90", "SBS91", "SBS92", "SBS93", "SBS94"))
  
  cesa = trinuc_mutation_rates(cesa, signature_exclusions = to_remove, signature_set = "COSMIC_v3.1",
                               signature_extractor = 'deconstructSigs')
  
  expect_error(trinuc_mutation_rates(cesa, signature_set = "COSMIC_v3.1"), 'have already been run')
  expect_error(trinuc_mutation_rates(cesa, signature_set = "COSMIC_v3.1", samples = cesa$samples[1]), 'have already been run')
  
  trinuc_ak = get_test_data("trinuc_mut_weighting.rds")
  expect_equal(length(waldo::compare(cesa@trinucleotide_mutation_weights, trinuc_ak, 
                                     tolerance = 1e-4, ignore_attr = 'index')), 0)
  
  
  # Ensure SBS counts (total and used by dS) look right
  expect_identical(cesa@trinucleotide_mutation_weights$signature_weight_table[, unlist(.(total_sbs, sig_extraction_sbs))],
                   c(67, 181, 38, 18, 0, 66, 180, 38, 0, 0))
  
  
  full_weight_table = get_signature_weights(cesa)
  expect_equal(full_weight_table[, .N], 5)
  expect_equal(full_weight_table[patient_id == "one_indel", unlist(.(total_sbs, sig_extraction_sbs))], c(0, 0))
  expect_equal(full_weight_table[patient_id == "zeroed-out", unlist(.(total_sbs, sig_extraction_sbs))], c(18, 0))
  
  # All tumors should have an above-threshold number of usable SBS, or they should be group-average-blended (but never both)
  expect_true(all(full_weight_table[, xor(sig_extraction_sbs > 49, group_avg_blended == T)]))
  
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
  raw_weight_input = prev_raw_sig[, .SD, .SDcols = patterns("(SBS)|(patient_id)")]
  cesa3 = set_signature_weights(cesa, signature_set = "COSMIC_v3.1", weights = raw_weight_input)
  
  
  which_pure = prev_relative_bio_sig[group_avg_blended == F, which = T]
  all.equal(cesa3@trinucleotide_mutation_weights$signature_weight_table[which_pure, -c("sig_extraction_sbs", "group_avg_blended")],
            prev_relative_bio_sig[which_pure, -c("sig_extraction_sbs", "group_avg_blended")])
  all.equal(cesa3@trinucleotide_mutation_weights$raw_signature_weights[,  .SD, .SDcols = patterns("(SBS)|(patient_id)")],
            raw_weight_input)
  all.equal(cesa3@trinucleotide_mutation_weights$trinuc_proportion_matrix[which_pure,],
            prev_rates_matrix[which_pure, ])
  
  # setting with the previous biological weights will yield identical weights/rates to the previous run
  cesa4 = set_signature_weights(cesa, signature_set = "COSMIC_v3.1", weights = prev_relative_bio_sig)
  expect_equal(cesa4@trinucleotide_mutation_weights$signature_weight_table[, -c("sig_extraction_sbs", "group_avg_blended")],
            prev_relative_bio_sig[, -c("sig_extraction_sbs", "group_avg_blended")])
  
  current = cesa4@trinucleotide_mutation_weights$trinuc_proportion_matrix
  expect_equal(current[sort(rownames(current)),], 
               prev_rates_matrix[sort(rownames(prev_rates_matrix)),])
  
  
  # Quick tests with MutationalPatterns (may be extended later)
  cesa5 = trinuc_mutation_rates(cesa, signatures_to_remove = to_remove, signature_set = "COSMIC_v3.1",
                                signature_extractor = "MutationalPatterns")
  
  # Compare ak rates (based off of deconstructSigs) with MP-derived rates
  mut_mat_1 = t(as.matrix(trinuc_ak$trinuc_proportion_matrix, rownames = 'patient_id'))
  mut_mat_2 = t(as.matrix(cesa5$trinuc_rates, rownames = 'patient_id'))
  expect_identical(colnames(mut_mat_1), colnames(mut_mat_2)) # ensure samples match up
  colnames(mut_mat_1) = paste0('dS_', colnames(mut_mat_1))
  colnames(mut_mat_2) = paste0('MP_', colnames(mut_mat_2))
  
  sim = MutationalPatterns::cos_sim_matrix(mut_mat_1, mut_mat_2)
  
  # The similarity matrix may be worth manually inspecting occasionally.
  # Here, we just verify the basic expectaton that each sample should have high similarity to itself,
  # except the dS_zeroed-out sample, which is not zeroed out in MP
  expect_gt(min(diag(sim)[1:4]) - .9, 0)
  expect_lt(sim[5,5], .7)
  
  # Ensure SBS counts (total and used by dS) look right
  # Difference from dS test above is that zeroed-out sample is not zeroed out in MP, so 18 usable SBS
  expect_identical(cesa5@trinucleotide_mutation_weights$signature_weight_table[, unlist(.(total_sbs, sig_extraction_sbs))],
                   c(67, 181, 38, 18, 0, 66, 180, 38, 18, 0))
  
  full_weight_table = get_signature_weights(cesa5)
  expect_equal(full_weight_table[, .N], 5)
  expect_equal(full_weight_table[patient_id == "one_indel", unlist(.(total_sbs, sig_extraction_sbs))], c(0, 0))
  # Note that it is c(18, 18) with use of MutationalPatters because of different calculation method
  expect_equal(full_weight_table[patient_id == "zeroed-out", unlist(.(total_sbs, sig_extraction_sbs))], c(18, 18))
  
  # All tumors should have an above-threshold number of usable SBS, or they should be group-average-blended (but never both)
  expect_true(all(full_weight_table[, xor(sig_extraction_sbs > 49, group_avg_blended == T)]))
})

test_that('trinuc_sbs_counts output is usable by extractors', {
  ds_counts = trinuc_sbs_counts(cesa$maf, genome = cesa$reference_data$genome, style = 'deconstructSigs')
  expect_silent(deconstructSigs::whichSignatures(tumor.ref = ds_counts, sample.id = 'good_A', contexts.needed = T, 
                                   signatures.ref = ces.refset.hg19$signatures$COSMIC_v3$signatures))
  mp_counts = trinuc_sbs_counts(cesa$maf, genome = cesa$reference_data$genome, style = 'MutationalPatterns')
  expect_silent(MutationalPatterns::fit_to_signatures(mut_matrix = mp_counts, signatures = t(ces.refset.hg19$signatures$COSMIC_v3$signatures)))
})


