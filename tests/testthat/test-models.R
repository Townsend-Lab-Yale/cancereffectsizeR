# will use pre-annotated CESAnalysis
cesa = load_cesa(get_test_file("cesa_hg38_for_test.rds"))


# Make sure baseline_mutation_rates is working (requires full mutation ID)
# Try taking just one AAC rate in one sample
rate_a = baseline_mutation_rates(cesa, variant_ids = "SCNN1D_G467R_ENSP00000368411.5", samples = "TCGA-BH-A18H")
rate_b = baseline_mutation_rates(cesa, aac_ids = "SCNN1D_G467R_ENSP00000368411.5", samples = "TCGA-BH-A18H")
expect_equal(rate_a, rate_b)
expect_equal(rate_a[[1]], 'TCGA-BH-A18H')
expect_equal(rate_a[[2]], 1.536797e-06)


# Try taking one SNV rate in one sample
snv_rate = baseline_mutation_rates(cesa, snv_ids = "1:151292923_C>T", 
                                     samples = 'P-0000015')
expect_equal(snv_rate[[1]], 'P-0000015')
expect_equal(snv_rate[[2]], 1.824055e-06)
snv_rate_other_way = baseline_mutation_rates(cesa, variant_ids = "1:151292923_C>T", 
                                             samples = 'P-0000015')
expect_identical(snv_rate, snv_rate_other_way)

# Invalid inputs
expect_error(baseline_mutation_rates(cesa), "no variants were input")
expect_error(baseline_mutation_rates(cesa, variant_ids = c("1:151292923_C>T", 'hi')), 
             "1 variant IDs are either invalid or not present")

# genes to plug into ces_variant; some with high-effect-size SNVs, others random
test_genes = c("EGFR", "ASXL3", "KRAS", "RYR2", "USH2A", "CSMD3", "TP53", 
               "CSMD1", "LRP1B", "ZFHX4", "FAT3", "CNTNAP5", "PCDH15", "NEB", 
               "RYR3", "PIK3CA", "ESR1", 
               "MAP3K1", "AKT1", "ARID1A", "FOXA1", 
               "TBX3", "PTEN")
main_effects_ak = fread(get_test_file("default_model_effects_brca_hg38.txt"))

test_that("ces_variant default effects", {
  cesa = ces_variant(cesa, variants = select_variants(cesa, genes = test_genes), run_name = 'main_test', cores = 1)
  main_effects = copy(cesa@selection_results[[1]])
  expect_equal(attr(main_effects, "si_cols"), "selection_intensity")
  expect_equal(main_effects[order(variant_id)], main_effects_ak[order(variant_id)], check.attributes = F, tolerance = 1e-4)
})


test_that("ces_variant with user-supplied variants", {
  # Test an SNV only covered in targeted panel.
  expected_si = ces_variant(cesa, variants = select_variants(cesa, variant_ids = '16:68823629_A>C'))@selection_results[[1]]$selection_intensity
  expect_equal(expected_si, 4187.333, tolerance = 1e-5)
  
  
  # Test an intergenic SNV. Should get same result regardless of hold-out option.
  expected_si = 7751.457
  expect_equal(ces_variant(cesa, variants = select_variants(cesa, variant_ids = "1:1468632_C>A"))@selection_results[[1]]$selection_intensity,
               expected_si, tolerance = 1e-4)
  expect_equal(ces_variant(cesa, variants = select_variants(cesa, variant_ids = "1:1468632_C>A"),
                           hold_out_same_gene_samples = F)@selection_results[[1]]$selection_intensity, 
               expected_si, tolerance = 1e-4)
})


test_that("ces_variant on subsets of samples", {
  # This variant appears appears 3 times in cherry, 1 in mountain apple, 0 in marionberry
  variant = select_variants(cesa, variant_ids = "FOXA1_F266L")
  
  just_cherry = ces_variant(cesa, variants = variant, samples = cesa$samples[fruit == 'cherry'])@selection_results[[1]]$selection_intensity
  also_marionberry = ces_variant(cesa, variants = variant, samples = cesa$samples[fruit %in% c('cherry', 'marionberry')])@selection_results[[1]]$selection_intensity
  with_all = ces_variant(cesa, variants = variant)@selection_results[[1]]$selection_intensity
  expect_lt(also_marionberry, just_cherry)
  expect_lt(also_marionberry, with_all)
  
  # Add some synthetic variants with no coverage; shouldn't affect selection
  new_maf = copy(cesa@maf)
  new_samples = copy(cesa@samples)
  new_maf[, Unique_Patient_Identifier := paste0(Unique_Patient_Identifier, '.1')] # cast as new samples
  new_samples[, covered_regions := "no_cov"]
  new_samples[, Unique_Patient_Identifier := paste0(Unique_Patient_Identifier, '.1')]
  new_rates = cesa@trinucleotide_mutation_weights$trinuc_proportion_matrix
  rownames(new_rates) = paste0(rownames(cesa@trinucleotide_mutation_weights$trinuc_proportion_matrix), '.1')
  without_cov = GenomicRanges::setdiff(cesa@coverage$exome$`exome+`, GRanges("14:37591985-37591990"))
  
  old_cesa = copy_cesa(cesa)
  cesa = add_covered_regions(cesa, covered_regions = without_cov, covered_regions_name = "no_cov",
                             coverage_type = "exome")
  aac_id = cesa@mutations$amino_acid_change[aac_id %like% 'FOXA1_F266L', aac_id]
  expect_identical(cesa@mutations$variants_to_cov[[aac_id]], c('exome+', 'top_genes'))
  
  # illegally adding new samples/data
  cesa@samples = rbind(cesa@samples, new_samples)
  cesa@maf = rbind(cesa@maf, new_maf[top_consequence != 'FOXA1_F266L']) # avoid adding uncovered records
  cesa@trinucleotide_mutation_weights$trinuc_proportion_matrix = rbind(cesa@trinucleotide_mutation_weights$trinuc_proportion_matrix,
                                                                       new_rates)
  
  # selection intensity should be unchanged by addition of new samples
  variant = select_variants(cesa, variant_ids = "FOXA1_F266L")
  all_with_new = ces_variant(cesa, variants = variant)@selection_results[[1]]
  expect_identical(with_all, all_with_new$selection_intensity)
  
  # this time, add data with ad-hoc genome-wide coverage; SI should go down
  cesa = add_covered_regions(cesa, covered_regions = cesa@coverage$exome$exome, covered_regions_padding = 1e9, 
                            coverage_type = "genome", covered_regions_name = 'whole_genome')
  cesa@samples[covered_regions == "no_cov", covered_regions := "whole_genome"]
  
  # Need to re-call since covered_in changed without samples changing. In the future,
  # the variants table should probably always be pulled from 
  #variant = select_variants(cesa, variant_ids = "EGFR_L858R_ENSP00000275493")
  with_wgs = ces_variant(cesa, variants = variant)@selection_results[[1]]
  expect_lt(with_wgs$selection_intensity, with_all)
  
  # Restore CESAnalysis for continued testing
  cesa = copy_cesa(old_cesa)
})

test_that("Variant-level epistasis", {
  results = ces_epistasis(cesa, variants = list(c("KRAS_G12V_ENSP00000256078.5", "GATA3 M293K")), conf = .9)@epistasis[[1]]
  to_test = results[, as.numeric(.(ces_A0, ces_B0, ces_A_on_B, ces_B_on_A, nA0,
                                   nB0, nAB, n00))]
  expect_equal(to_test[1:4], c(13057.7094428118, 11974.1862697053, 0.001, 0.001), tolerance = 1e-5)
  expect_equal(to_test[5:8], c(7, 6, 0, 1077))
  ci = as.numeric(results[, .SD, .SDcols = patterns("ci")])
  expect_equal(ci, c(6527.46003622592, 22938.2519029567, 5623.19405835833, 21906.9674813457, 
                     NA, 897584.501451215, NA, 744152.729482514), tolerance = 1e-5)
})


test_that("Gene-level SNV epistasis analysis", {
  variants_to_use = cesa$variants[gene %in% c('EGFR', 'KRAS', 'TP53') & samples_covering == cesa$samples[, .N]]
  cesa = ces_gene_epistasis(cesa, genes = c("EGFR", "KRAS", "TP53"), variants = variants_to_use, conf = .95)
  results_ak = get_test_data("epistatic_effects.rds")
  expect_equal(cesa@epistasis[[1]], results_ak, tolerance = 1e-5)
  
  # now simulate with compound variants, should get almost same results
  comp = define_compound_variants(cesa, 
                                  variant_table = select_variants(cesa, genes = c("EGFR", "KRAS", "TP53"))[samples_covering == cesa$samples[, .N]],
                                  by = "gene", merge_distance = Inf)
  cesa = ces_epistasis(cesa, comp, conf = .95)
  all.equal(cesa@epistasis[[1]][, -c(1, 2)], cesa@epistasis[[2]][, -c(1, 2)], check.attributes = F, tolerance = 1e-4)
  
  # variant counts should always add up
  expect_equal(cesa$epistasis[[1]][, unique(nA0 + nB0 + nAB + n00 - n_total)], 0)
  
  # For the purposes of this test, we need a gene with variable coverage and one just covered in WXS.
  # These first two tests are to make sure these condiditons aren't accidentally violated if test data changes in the future.
  expect_gt(cesa$variants[gene == 'ARID1A', uniqueN(samples_covering)], 1)
  expect_equal(cesa$variants[gene == 'TTN', unique(samples_covering)], cesa$samples[covered_regions == 'exome+', .N])
  cesa = expect_message(expect_warning(ces_gene_epistasis(cesa = cesa, genes = c('ARID1A', 'TTN'), variants = cesa$variants[gene %in% c('ARID1A', 'TTN')], 
                            run_name = 'early_output'), 'this variant pair had to be skipped'), 'all NAs')
  early_output = cesa$epistasis$early_output
  expect_equal(early_output[, unique(c(nA0, nB0, nAB, n00, n_total))], 0)
  
  # 4 parameters and 8 low/high CIs should all be NA
  expect_equal(early_output[, as.numeric(.SD), .SDcols = patterns('ces')], rep(NA_real_, 12))
  
  covered_variants = select_variants(cesa = cesa, genes = c('ARID1A', 'TTN'), gr = cesa$coverage_ranges$exome$`exome+`)
  cesa = ces_gene_epistasis(cesa, genes = c('ARID1A', 'TTN'), variants = covered_variants, run_name = 'should_work')
  expect_true(cesa$epistasis$should_work[, ! anyNA(.(ces_A0, ces_B0, ces_A_on_B, ces_B_on_A))])
  
  # Should get a note that no eligible variants have joint coverage. Unlike the above check of ARID1A/TTN,
  # there won't be a warning about skipping the variant.
  expect_message(ces_gene_epistasis(cesa = cesa, genes = c('KRAS', 'ARID1A'), samples = cesa$samples[1:100]),
                 'are all NAs')
  
})

test_that("Compound variant creation", {
  recurrents = select_variants(cesa, min_freq = 2)
  kras_12_13 = recurrents[gene == "KRAS" & aa_pos %in% c(12, 13)]
  expect_equal(kras_12_13[, .N], 4)
  all_kras_12_13 = select_variants(cesa, variant_position_table = kras_12_13)
  expect_equal(all_kras_12_13[, .N], 6)
  
  # Should recapitulate single variant results with 1-SNV-size compound variants
  # Note, though, that compound variant mutation rates get calculated via the SNV method, which 
  # averages any overlapping gene rates, so sometimes rates will vary between a one-SNV AAC and the equivalent compound variant.
  single_snv_comp = define_compound_variants(cesa, all_kras_12_13, by = c("gene", "aa_alt"), merge_distance = 0)
  expect_equal(length(single_snv_comp), 6)
  results = ces_variant(cesa, variants = single_snv_comp, hold_out_same_gene_samples = T)
  same_results = ces_variant(cesa, variants = single_snv_comp, hold_out_same_gene_samples = T)
  expect_equal(same_results@selection_results, results@selection_results)

  # this time, do per position
  pos_comp = define_compound_variants(cesa, all_kras_12_13, merge_distance = 0)
  expect_equal(length(pos_comp), 2)
  expect_equal(pos_comp$snv_info[, .N], 6)
  results = ces_variant(cesa, variants = pos_comp, hold_out_same_gene_samples = T, conf = NULL)$selection[[1]]
  expect_equal(results$selection_intensity, c(3545.80843024223, 5950.80593547194), tolerance = 1e-5)
  
  expect_equal(length(define_compound_variants(cesa, recurrents, merge_distance = 1000)), 52)
  expect_equal(length(define_compound_variants(cesa, recurrents, merge_distance = 1e9)), 15)
  expect_equal(length(suppressWarnings(define_compound_variants(cesa, recurrents, merge_distance = Inf))), 1)
  expect_equal(length(define_compound_variants(cesa, recurrents, by = "aa_alt", merge_distance = Inf)), 22)
  
  # causes some but not all same-gene mutations to get merged
  expect_equal(length(define_compound_variants(cesa, recurrents, by = c("gene","aa_ref"), merge_distance = 1000)), 97)
})

test_that("Attribution to mutational sources", {
  mut_effects = mutational_signature_effects(cesa = cesa, effects = main_effects_ak, 
                                             samples = cesa$samples[coverage == 'exome'])
  mut_ak = get_test_data('mut_effect_attribution.rds')
  expect_equal(mut_effects, mut_ak)
})


test_that("Just genome data", {
  # Various edge cases being hit: Just WGS data, load_maf() on only previously annotated variants,
  # running ces_variant() on just 2 samples, running with SNVs and indels but no DBS.
  c2 = CESAnalysis(ces.refset.hg38)
  c2 = add_variants(target_cesa = c2, source_cesa = cesa)
  c2 = load_maf(c2, maf = cesa$maf[1:100], coverage = 'genome')
  c2 = assign_group_average_trinuc_rates(c2)
  c2 = set_gene_rates(c2, rates = cesa$gene_rates[, .(pid, rate = rate_grp_1)])
  c2 = ces_variant(c2, variants = c2$variants[c2$maf$variant_id[1:5], on = 'variant_id'])
  expect_equivalent(unique(c2$selection$selection.1$included_with_variant), 1)
  expect_equivalent(unique(c2$selection$selection.1$included_total), 2)
})





