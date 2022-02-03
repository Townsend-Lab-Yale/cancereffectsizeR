# will use pre-annotated CESAnalysis
cesa = load_cesa(get_test_file("annotated_fruit_cesa.rds"))
trinuc_rates = fread(get_test_file("luad_hg19_trinuc_rates.txt"))
cesa = set_trinuc_rates(cesa, trinuc_rates = trinuc_rates)
sig_weights = fread(get_test_file("luad_hg19_sig_table_biological.txt"))
cesa@trinucleotide_mutation_weights$signature_weight_table = sig_weights
cesa@advanced$snv_signatures[["COSMIC v3.1"]] = get_ces_signature_set("ces.refset.hg19", "COSMIC_v3.1")

precalc_rates = fread(get_test_file("luad_fruit_gene_rates.txt"))
cesa = set_gene_rates(cesa, precalc_rates[, .(gene, rate_grp_1)], samples = cesa$samples[group == 'marionberry'])
cesa = set_gene_rates(cesa, precalc_rates[, .(gene, rate_grp_2)], samples = cesa$samples[c("cherry", "mountain_apple"), on = 'group'])

# Make sure baseline_mutation_rates is working (requires full mutation ID)
# Try taking just one AAC rate in one sample
# Try taking one AAC rate
rate_a = baseline_mutation_rates(cesa, variant_ids = "FAT1_D1508H_ENSP00000406229", samples = "sample-5")
rate_b = baseline_mutation_rates(cesa, aac_ids = "FAT1_D1508H_ENSP00000406229", samples = "sample-5")
expect_equal(rate_a, rate_b)
expect_equal(rate_a[[1]], 'sample-5')
expect_equal(rate_a[[2]], 2.843945e-06)


# Try taking one SNV rate in one sample
snv_rate = baseline_mutation_rates(cesa, snv_ids = "3:186276363_T>A", 
                                     samples = 'sample-106')
expect_equal(snv_rate[[1]], 'sample-106')
expect_equal(snv_rate[[2]], 1.20681e-06)
snv_rate_other_way = baseline_mutation_rates(cesa, variant_ids = "3:186276363_T>A", 
                                             samples = 'sample-106')
expect_identical(snv_rate, snv_rate_other_way)

# Invalid inputs
expect_error(baseline_mutation_rates(cesa), "no variants were input")
expect_error(baseline_mutation_rates(cesa, variant_ids = c("3:186276363_T>A", 'hi')), 
             "1 variant IDs are either invalid or not present")

# genes to plug into ces_variant; some with high-effect-size SNVs, others random
test_genes = c("EGFR", "ASXL3", "KRAS", "RYR2", "USH2A", "CSMD3", "TP53", "CSMD1", "LRP1B",
               "ZFHX4", "FAT3", "CNTNAP5", "PCDH15", "NEB", "RYR3", "DMD", "KATNAL1",
               "OR13H1", "KSR1")
test_that("ces_variant with sswm", {
  cesa = ces_variant(cesa, variants = select_variants(cesa, genes = test_genes), cores = 1)
  results = cesa@selection_results[[1]]
  expect_equal(attr(results, "si_cols"), "selection_intensity")
  results_ak = fread(get_test_file("fruit_sswm_out.txt"))
  expect_equal(results[order(variant_id)], results_ak[order(variant_id)], check.attributes = F, tolerance = 1e-7)
})

test_that("ces_variant with sswm_sequential", {
  cesa = expect_warning(ces_variant(cesa, variants = select_variants(cesa, genes = c("EGFR", "KRAS", "TP53"), variant_ids = "CR2 R247L"), 
                                    model = "sswm_sequential", groups = list(c("marionberry", "cherry"), "mountain_apple")),
                        'groups is deprecated')
  results = cesa@selection_results[[1]] # previous run not saved due to test_that scoping
  expect_equal(attr(results, "si_cols"), c("si_1", "si_2"))
  results_ak = fread(get_test_file("fruit_sswm_sequential_out.txt"))
  expect_equal(results[order(variant_id)], results_ak[order(variant_id)], check.attributes = F, tolerance = 1e-6)
})

test_that("ces_variant bad groups inputs", {
  expect_error(suppressWarnings(ces_variant(cesa, model = "sswm_sequential", groups = c("cherry","cherry"))), "groups are re-used")
  expect_error(ces_variant(cesa, model = "sequential", ordering_col = 'group', ordering = c("cherry","cherry")),
               "ordering contains repeated values")
  expect_error(suppressWarnings(ces_variant(cesa, model = "sswm_sequential", groups = list(c("cherry", "marionberry")))),
               "should be a list with length at least two")
  expect_message(ces_variant(cesa, variants = select_variants(cesa, variant_ids = "CR2 R247L"),
                             samples = cesa$samples[1:10]),
                  "samples are being excluded from selection inference")
})

test_that("ces_variant with user-supplied variants", {
  expect_equal(ces_variant(cesa, variants = select_variants(cesa, variant_ids = "10:100190376_C>A"))@selection_results[[1]][, .N], 1)
  variants = select_variants(cesa, genes = c("TP53", "KRAS"), variant_ids = "KRAS_G12D_ENSP00000256078", min_freq = 2)
  expect_equal(ces_variant(cesa, variants = variants)@selection_results[[1]][, .N], 9)
  
  # Test an SNV that doesn't overlap a gene. Should get same result regardless of hold-out option.
  expected_si = 5134.611
  expect_equal(ces_variant(cesa, variants = select_variants(cesa, variant_ids = '1:1903518_G>A'))@selection_results[[1]]$selection_intensity,
               expected_si, tolerance = 1e-4)
  expect_equal(ces_variant(cesa, variants = select_variants(cesa, variant_ids = '1:1903518_G>A'),
                           hold_out_same_gene_samples = F)@selection_results[[1]]$selection_intensity, 
               expected_si, tolerance = 1e-4)
})


#baseline_mutation_rates(cesa, aac_ids = 'EGFR_L858R_ENSP00000275493')
test_that("ces_variant on subsets of samples", {
  # EGFR L858R appears 5 times in cherry, 4 in mountain apple, 0 in marionberry
  egfr = select_variants(cesa, variant_ids = "EGFR_L858R_ENSP00000275493")
  
  just_cherry = ces_variant(cesa, variants = egfr, samples = cesa$samples[group == 'cherry'])@selection_results[[1]]$selection_intensity
  also_marionberry = expect_warning(ces_variant(cesa, variants = egfr, 
                                                groups = c("cherry", "marionberry"))@selection_results[[1]]$selection_intensity,
                                    'groups is deprecated')
  with_all = ces_variant(cesa, variants = egfr)@selection_results[[1]]$selection_intensity
  expect_lt(also_marionberry, just_cherry)
  expect_lt(also_marionberry, with_all)
  
  
  # add some "new" data with no coverage at EGFR L858R
  # Add some synthetic variants with no coverage; shouldn't affect selection
  new_maf = copy(cesa@maf)
  new_samples = copy(cesa@samples)
  new_maf[, Unique_Patient_Identifier := paste0(Unique_Patient_Identifier, '.1')] # cast as new samples
  new_samples[, covered_regions := "no_egfr"]
  new_samples[, Unique_Patient_Identifier := paste0(Unique_Patient_Identifier, '.1')]
  new_rates = cesa@trinucleotide_mutation_weights$trinuc_proportion_matrix
  rownames(new_rates) = paste0(rownames(cesa@trinucleotide_mutation_weights$trinuc_proportion_matrix), '.1')
  no_egfr_cov = GenomicRanges::setdiff(cesa@coverage$exome$`exome+`, GRanges("7:55259514-55259516"))
  
  cesa = add_covered_regions(cesa, covered_regions = no_egfr_cov, covered_regions_name = "no_egfr",
                             coverage_type = "exome")
  
  egfr_variants = select_variants(cesa, genes = "EGFR", min_freq = 2)
  expect_true(egfr_variants[variant_name == "EGFR_L858R", identical(unlist(covered_in), "exome+")])
  expect_true(egfr_variants[variant_name != "EGFR_L858R", length(unlist(unique(covered_in))) == 2])
  
  # illegally adding new samples/data
  cesa@samples = rbind(cesa@samples, new_samples)
  cesa@maf = rbind(cesa@maf, new_maf[variant_id != '7:55259515_T>G']) # avoid adding uncovered record
  cesa@trinucleotide_mutation_weights$trinuc_proportion_matrix = rbind(cesa@trinucleotide_mutation_weights$trinuc_proportion_matrix,
                                                                       new_rates)
  
  # selection intensity should be unchanged by addition of new samples
  all_with_new = ces_variant(cesa, variants = egfr)@selection_results[[1]]
  expect_identical(with_all, all_with_new$selection_intensity)
  
  # this time, add data with ad-hoc genome-wide coverage; SI should go down
  cesa = add_covered_regions(cesa, covered_regions = cesa@coverage$exome$exome, covered_regions_padding = 1e9, 
                            coverage_type = "genome", covered_regions_name = 'whole_genome')
  cesa@samples[covered_regions == "no_egfr", covered_regions := "whole_genome"]
  
  # need to re-call since covered_in changed without samples changing
  egfr = select_variants(cesa, variant_ids = "EGFR_L858R_ENSP00000275493")
  with_wgs = ces_variant(cesa, variants = egfr)@selection_results[[1]]
  expect_lt(with_wgs$selection_intensity, with_all)
})

test_that("Variant-level epistasis", {
  results = ces_epistasis(cesa, variants = list(c("KRAS G12V", "GSTP1_L184L")), conf = .9)@epistasis[[1]]
  to_test = results[, as.numeric(.(ces_v1, ces_v2, ces_v1_after_v2, ces_v2_after_v1, joint_cov_samples_just_v1,
                                   joint_cov_samples_just_v2, joint_cov_samples_with_both, joint_cov_samples_with_neither))]
  expect_equal(to_test[1:4], c(26867.5225688425, 1441.71206944064, 281837.785451727, 13805.1800777728), tolerance = 1e-3)
  expect_equal(to_test[5:8], c(6, 2, 1, 100))
  ci = as.numeric(results[, .SD, .SDcols = patterns("ci")])
  expect_equal(ci, c(12795.1211307558, 48488.3936128441, 362.207984931453, 3610.84903759668, 
                     NA, 2180128.54763536, NA, 88952.7367254685), tolerance = 1e-3)
})


test_that("Gene-level SNV epistasis analysis", {
  cesa = ces_gene_epistasis(cesa, genes = c("EGFR", "KRAS", "TP53"), variants = 'recurrent', conf = .95)
  results_ak = get_test_data("epistasis_results.rds")
  expect_equal(cesa@epistasis[[1]], results_ak, tolerance = 1e-3)
  
  # now simulate with compound variants, should get almost same results
  comp = define_compound_variants(cesa, variant_table = select_variants(cesa, genes = c("EGFR", "KRAS", "TP53"), min_freq = 2),
                                  by = "gene", merge_distance = Inf)
  cesa = ces_epistasis(cesa, comp, conf = .95)
  all.equal(cesa@epistasis[[1]][, -c(1, 2)], cesa@epistasis[[2]][, -c(1, 2)], check.attributes = F, tolerance = 1e-4)
})

test_that("Compound variant creation", {
  recurrents = select_variants(cesa, min_freq = 2)
  kras_12_13 = recurrents[gene == "KRAS" & aa_pos %in% c(12, 13)]
  expect_equal(kras_12_13[, .N], 4)
  all_kras_12_13 = select_variants(cesa, variant_position_table = kras_12_13)
  expect_equal(all_kras_12_13[, .N], 5)  # one has MAF freq of 1
  
  # should recapitulate single variant results with 1-SNV-size compound variants
  single_snv_comp = define_compound_variants(cesa, all_kras_12_13, by = c("gene", "aa_alt"), merge_distance = 0)
  expect_equal(length(single_snv_comp), 5)
  results = expect_warning(ces_variant(cesa, variants = single_snv_comp, model = "sswm_sequential", 
                                         groups = list(c("marionberry", "cherry"), "mountain_apple"), hold_out_same_gene_samples = T),
                             'groups is deprecated')

  same_results = ces_variant(cesa, variants = single_snv_comp, model = 'sequential', ordering_col = 'group',
                             ordering = list(c("marionberry", "cherry"), "mountain_apple"), hold_out_same_gene_samples = T)
  expect_equal(same_results@selection_results, results@selection_results)
  results = results$selection[[1]][c("KRAS.Cys.1", "KRAS.Asp.1", "KRAS.Cys.2", "KRAS.Asp.2", "KRAS.Val.1"), on = "variant_name"]
  prev = fread(get_test_file("fruit_sswm_sequential_out.txt"))[all_kras_12_13$variant_id, on = "variant_id"][, 3:4]
  expect_equal(results[, 3:4], prev, tolerance = 1e-6)
  
  # this time, do per position
  pos_comp = define_compound_variants(cesa, all_kras_12_13, merge_distance = 0)
  expect_equal(length(pos_comp), 2)
  expect_equal(pos_comp$snv_info[, .N], 5)
  results = ces_variant(cesa, variants = pos_comp, hold_out_same_gene_samples = T, conf = NULL)$selection[[1]]
  expect_equal(results$selection_intensity, c(6482.0086263444, 29946.7385123648), tolerance = 1e-3)
  # If above output changes, be sure SIs remain between lowest/highest single-variant SIs at each site (see .rds)
  
  expect_equal(length(define_compound_variants(cesa, recurrents, merge_distance = 1000)), 40)
  expect_equal(length(define_compound_variants(cesa, recurrents, merge_distance = 1e9)), 19)
  expect_equal(length(define_compound_variants(cesa, recurrents, merge_distance = Inf)), 1)
  expect_equal(length(define_compound_variants(cesa, recurrents, by = "aa_alt", merge_distance = Inf)), 19)
  
  # causes two TP53 Arg mutations to get merged
  expect_equal(length(define_compound_variants(cesa, recurrents, by = c("gene","aa_ref"), merge_distance = 1000)), 43)
})








