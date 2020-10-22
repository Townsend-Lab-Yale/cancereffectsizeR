
test_that("gene mutation rates", {
  dndscv_cesa = load_cesa(get_test_file("annotated_fruit_cesa.rds"))
  dndscv_input = get_test_data("dndscv_input_single.rds")
  dndscv_raw_output = get_test_data("dndscv_raw_output_single.rds")

  # Run gene_mutation_rates, but skip running dNdScv. Instead, verify that input to
  # run_dndscv matches what's expected, then load pre-generated output and continue
  dndscv_cesa = mockr::with_mock(
    run_dndscv = function(mutations, gene_list, cv, gr_genes, refdb) {
      expect_equal(list(mutations = mutations, gene_list = gene_list, cv = cv), dndscv_input)
      return(dndscv_raw_output)
    },
    gene_mutation_rates(dndscv_cesa, covariates = "lung"))
  rates_ak = get_test_data("mutrates.rds")
  expect_equal(dndscv_cesa@mutrates, rates_ak)
})

# will use pre-annotated CESAnalysis
cesa = load_cesa(get_test_file("annotated_fruit_cesa.rds"))
trinuc_rates = fread(get_test_file("luad_hg19_trinuc_rates.txt"))
cesa = set_trinuc_rates(cesa, trinuc_rates = trinuc_rates)
sig_weights = fread(get_test_file("luad_hg19_sig_table.txt"))
cesa@trinucleotide_mutation_weights$signature_weight_table = sig_weights
cesa@advanced$snv_signatures = get_ces_signature_set("ces_hg19_v1", "COSMIC_v3.1")

precalc_rates = fread(get_test_file("luad_fruit_gene_rates.txt"))
cesa = set_gene_rates(cesa, precalc_rates[, .(gene, rate_grp_1)], sample_group = "marionberry")
cesa = set_gene_rates(cesa, precalc_rates[, .(gene, rate_grp_2)], sample_group = c("cherry", "mountain_apple"))


# genes to plug into ces_variant; some with high-effect-size SNVs, others random
test_genes = c("EGFR", "ASXL3", "KRAS", "RYR2", "USH2A", "CSMD3", "TP53", "CSMD1", "LRP1B",
               "ZFHX4", "FAT3", "CNTNAP5", "PCDH15", "NEB", "RYR3", "DMD", "KATNAL1",
               "OR13H1", "KSR1")
test_that("ces_variant with sswm", {
  cesa = ces_variant(cesa, variants = select_variants(cesa, genes = test_genes), cores = 1)
  results = cesa@selection_results[[1]]
  results_ak = fread(get_test_file("fruit_sswm_out.txt"))
  expect_equal(results[order(variant_id)], results_ak[order(variant_id)])
})

test_that("ces_variant with sswm_sequential", {
  cesa = ces_variant(cesa, variants = select_variants(cesa, genes = c("EGFR", "KRAS", "TP53"), variant_passlist = "CR2 R247L"),
                 lik_fn = "sswm_sequential", group_ordering = list(c("marionberry", "cherry"), "mountain_apple"))
  results = cesa@selection_results[[1]] # previous run not saved due to test_that scoping
  results_ak = fread(get_test_file("fruit_sswm_sequential_out.txt"))
  expect_equal(results[order(variant_id)], results_ak[order(variant_id)])
})

test_that("ces_variant bad group_ordering inputs", {
  expect_error(ces_variant(cesa, lik_fn = "sswm_sequential", group_ordering = c("cherry","cherry")), "groups are re-used")
  expect_error(ces_variant(cesa, lik_fn = "sswm_sequential", group_ordering = list(c("cherry", "marionberry"))),
               "must have length of at least 2")
  expect_warning(ces_variant(cesa, variants = select_variants(cesa, variant_passlist = "CR2 R247L"),
                             lik_fn = "sswm_sequential", group_ordering = c("cherry", "marionberry")),
                  "The following groups were not included in group_ordering")
})

test_that("ces_variant with user-supplied variants", {
  expect_equal(ces_variant(cesa, variants = select_variants(cesa, variant_passlist = "10:100190376_C>A"))@selection_results[[1]][, .N], 1)
  variants = select_variants(cesa, genes = c("TP53", "KRAS"), variant_passlist = "KRAS_G12D_ENSP00000256078", min_freq = 2)
  expect_equal(ces_variant(cesa, variants = variants)@selection_results[[1]][, .N], 9)
})

test_that("Gene-level SNV epistasis analysis", {
  results = ces_gene_epistasis(cesa, genes = c("EGFR", "KRAS", "TP53"), conf = .95)
  results_ak = get_test_data("epistasis_results.rds")
  expect_equal(results, results_ak, tolerance = 1e-3)
})

test_that("Variant-level epistasis", {
  results = ces_epistasis(cesa, variants = list(c("KRAS G12V", "GSTP1_L184L")), conf = .9)
  to_test = results[, as.numeric(.(ces_v1, ces_v2, ces_v1_after_v2, ces_v2_after_v1, covered_tumors_just_v1,
                                   covered_tumors_just_v2, covered_tumors_with_both, covered_tumors_with_neither))]
  expect_equal(to_test, c(28864.82, 1432.69, 549695.3, 9508.182, 6, 2, 1, 100), tolerance = 1e-3)
  ci = as.numeric(results[, .SD, .SDcols = patterns("ci")])
  expect_equal(ci, c(13629.84, 52477.08, 382.9144, 3470.92, NA, 3610458, NA, 85519.34), tolerance = 1e-3)
})






