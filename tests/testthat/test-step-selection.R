# Tests for step-specific selection: step_selection_lik(), assign_stage_index(),
# stage_mutation_proportions(), ces_variant_step(), step_selection_LRT(), plot_effects_step().

## ---- step_selection_lik() should collapse to sswm_lik() under degenerate proportions ----
test_that("step_selection_lik degenerates correctly to sswm_lik", {
  set.seed(1)
  rates_with = setNames(runif(5, 1e-7, 1e-6), paste0("S", 1:5))
  rates_without = setNames(runif(20, 1e-7, 1e-6), paste0("S", 6:25))
  all_samples = c(names(rates_with), names(rates_without))

  # Everyone in stage 2 (group_index = 2); with stage_mut_prop = c(0, 1), gamma[1] (stage 1's SI)
  # should have zero influence, and the likelihood should exactly match sswm_lik()'s single-SI model
  # evaluated at gamma[2].
  sample_index = rbind(
    data.table(Unique_Patient_Identifier = all_samples, group_index = 2L, group_name = "late"),
    data.table(Unique_Patient_Identifier = "placeholder_stage1_sample", group_index = 1L, group_name = "early")
  )

  fn_step = step_selection_lik(rates_with, rates_without, sample_index, stage_mut_prop = c(0, 1))
  fn_basic = sswm_lik(rates_with, rates_without)

  for (g2 in c(500, 1000, 5000, 50000)) {
    expect_equal(fn_step(c(g2 * 3.7, g2)), fn_basic(g2), tolerance = 1e-8)
  }

  expect_equal(bbmle::parnames(fn_step), c("si_early", "si_late"))
  expect_error(step_selection_lik(rates_with, rates_without, sample_index, stage_mut_prop = c(0, .5, .5)),
              "must match")
})

## ---- assign_stage_index() ----
test_that("assign_stage_index builds a correct sample_index", {
  cesa = load_cesa(get_test_file("cesa_hg38_for_test.rds"))

  si = assign_stage_index(cesa, stage_col = "fruit",
                          stage_order = list(early = "marionberry", mid = "cherry", late = "mountain_apple"))
  expect_equal(si[, .N], cesa$samples[, .N])
  expect_setequal(names(si), c("Unique_Patient_Identifier", "group_index", "group_name"))
  expect_equal(si[group_name == "early", uniqueN(group_index)], 1)
  expect_equal(si[group_name == "early", group_index[1]], 1L)
  expect_equal(si[group_name == "late", group_index[1]], 3L)

  # plain vector form: stage names default to the values themselves
  si2 = assign_stage_index(cesa, stage_col = "fruit", stage_order = c("marionberry", "cherry", "mountain_apple"))
  expect_setequal(si2$group_name, c("marionberry", "cherry", "mountain_apple"))

  expect_error(assign_stage_index(cesa, stage_col = "fruit", stage_order = c("marionberry", "cherry")),
              "not covered by stage_order")
  expect_error(assign_stage_index(cesa, stage_col = "not_a_column", stage_order = c("a", "b")),
              "not a column")
  expect_error(assign_stage_index(cesa, stage_col = "fruit", stage_order = "marionberry"),
              "at least 2 stages")
})

## ---- stage_mutation_proportions(), including the monotonicity guard ----
test_that("stage_mutation_proportions computes correct proportions and flags non-monotonic genes", {
  cesa = load_cesa(get_test_file("cesa_hg38_for_test.rds"))
  cesa = clear_gene_rates(cesa)

  genes = c("EGFR", "ASXL3", "KRAS", "RYR2", "TP53")
  rate1 = data.table(gene = genes, rate = c(1e-6, 2e-6, 1.5e-6, 3e-6, 5e-7))
  # RYR2's stage-2 cumulative rate (2e-6) is deliberately below stage 1's (3e-6): invalid/non-monotonic
  rate2 = data.table(gene = genes, rate = c(2e-6, 3e-6, 2.5e-6, 2e-6, 1.2e-6))

  all_samples = cesa$samples$Unique_Patient_Identifier
  half = ceiling(length(all_samples) / 2)
  grpA = all_samples[1:half]
  grpB = all_samples[(half + 1):length(all_samples)]

  cesa = set_gene_rates(cesa, rates = rate1, samples = grpA, missing_genes_take_nearest = TRUE)
  cesa = set_gene_rates(cesa, rates = rate2, samples = grpB, missing_genes_take_nearest = TRUE)

  expect_warning(
    smp <- stage_mutation_proportions(cesa, rate_cols = c("rate_grp_1", "rate_grp_2"),
                                      stage_names = c("early", "late"), rate_source = "gene_rates"),
    "not non-decreasing"
  )

  egfr_row = smp[gene == "EGFR"]
  expect_equal(egfr_row$p_early, 0.5, tolerance = 1e-9)
  expect_equal(egfr_row$p_late, 0.5, tolerance = 1e-9)
  expect_equal(egfr_row$p_early + egfr_row$p_late, 1, tolerance = 1e-9)

  ryr2_row = smp[gene == "RYR2"]
  expect_equal(ryr2_row$p_early + ryr2_row$p_late, 1, tolerance = 1e-9)
  expect_true(ryr2_row$p_late > 0) # floored to a small positive share, not zero
  expect_equal(egfr_row$rate_early, 1e-6)
  expect_equal(egfr_row$rate_late, 2e-6)

  expect_error(
    stage_mutation_proportions(cesa, rate_cols = c("rate_grp_1", "rate_grp_2"), rate_source = "gene_rates", on_invalid = "error"),
    "not non-decreasing"
  )

  smp_na = suppressWarnings(stage_mutation_proportions(cesa, rate_cols = c("rate_grp_1", "rate_grp_2"), rate_source = "gene_rates", on_invalid = "NA"))
  expect_true(all(is.na(smp_na[gene == "RYR2", .(p_rate_grp_1, p_rate_grp_2)])))
  expect_false(anyNA(smp_na[gene == "EGFR"]))

  expect_error(stage_mutation_proportions(cesa, rate_cols = c("rate_grp_1"), rate_source = "gene_rates"), "at least two")
  expect_error(stage_mutation_proportions(cesa, rate_cols = c("rate_grp_1", "not_a_col"), rate_source = "gene_rates"), "not present")

  # default rate_source needs full dNdScv output, which set_gene_rates() doesn't produce
  expect_error(stage_mutation_proportions(cesa, rate_cols = c("rate_grp_1", "rate_grp_2")), "not present|Full dNdScv output")
})

test_that("stage_mutation_proportions uses uncorrected dNdScv expected rates by default", {
  cesa = load_cesa(get_test_file("cesa_hg38_for_test.rds"))
  genes = c("EGFR", "KRAS", "TP53")
  # Minimal stand-ins for full dNdScv output: expected synonymous counts and mutated samples
  fake_dndscv = function(exp_syn, num_samples) {
    list(genemuts = data.frame(gene_name = genes, exp_syn_cv = exp_syn, n_syn = c(50, 50, 50)),
         annotmuts = data.frame(sampleID = paste0("s", seq_len(num_samples))))
  }
  cesa@dndscv_out_list = list(rate_grp_1 = fake_dndscv(c(1, 2, 3), 10),
                              rate_grp_2 = fake_dndscv(c(8, 8, 9), 20))
  smp = stage_mutation_proportions(cesa, stage_names = c("early", "late"))

  # synonymous site counts cancel in the proportions: p_early = (e1 / n1) / (e2 / n2);
  # observed synonymous counts (n_syn) must not matter
  expect_equal(smp[genes, p_early, on = "gene"], c(1/10, 2/10, 3/10) / (c(8, 8, 9) / 20))
  expect_equal(smp$p_early + smp$p_late, rep(1, 3))

  RefCDS = cancereffectsizeR:::.ces_ref_data[[cesa@ref_key]]$RefCDS.dndscv
  if (is.null(RefCDS)) RefCDS = cancereffectsizeR:::.ces_ref_data[[cesa@ref_key]]$RefCDS
  sites = sapply(RefCDS[genes], function(x) colSums(x[["L"]])[1])
  expect_equal(smp[genes, rate_late, on = "gene"], unname(c(8, 8, 9) / sites / 20))

  cesa@dndscv_out_list = list(rate_grp_1 = data.table(gene_name = genes), rate_grp_2 = data.table(gene_name = genes))
  expect_error(stage_mutation_proportions(cesa), "save_all_dndscv_output")
})

## ---- ces_variant_step() end to end, plus step_selection_LRT() and plot_effects_step() ----
test_that("ces_variant_step produces expected output and integrates with step_selection_LRT/plot_effects_step", {
  cesa = load_cesa(get_test_file("cesa_hg38_for_test.rds"))
  cesa = clear_gene_rates(cesa)

  genes = c("EGFR", "ASXL3", "KRAS", "TP53")
  rate1 = data.table(gene = genes, rate = c(1e-6, 2e-6, 1.5e-6, 3e-6))
  rate2 = data.table(gene = genes, rate = c(2e-6, 3e-6, 2.5e-6, 5e-6))

  all_samples = cesa$samples$Unique_Patient_Identifier
  half = ceiling(length(all_samples) / 2)
  grpA = all_samples[1:half]
  grpB = all_samples[(half + 1):length(all_samples)]

  cesa = set_gene_rates(cesa, rates = rate1, samples = grpA, missing_genes_take_nearest = TRUE)
  cesa = set_gene_rates(cesa, rates = rate2, samples = grpB, missing_genes_take_nearest = TRUE)

  sample_index = data.table(
    Unique_Patient_Identifier = c(grpA, grpB),
    group_index = c(rep(1L, length(grpA)), rep(2L, length(grpB))),
    group_name = c(rep("early", length(grpA)), rep("late", length(grpB)))
  )
  setkey(sample_index, "Unique_Patient_Identifier")

  # genome-wide nearest-neighbor imputation (missing_genes_take_nearest = TRUE, above) means most
  # genes besides our 4 targets get non-monotonic proportions purely from imputation noise; that's
  # expected here and covered by its own test above, so just suppress the warning for this one.
  smp = suppressWarnings(stage_mutation_proportions(cesa, rate_cols = c("rate_grp_1", "rate_grp_2"),
                                                    stage_names = c("early", "late"), rate_source = "gene_rates"))

  kras_variants = select_variants(cesa, genes = "KRAS", min_freq = 2)
  expect_true(kras_variants[, .N] > 0)

  # per-stage rates still in place: should warn
  expect_warning(ces_variant_step(cesa, variants = kras_variants, stage_mut_prop = smp, sample_index = sample_index,
                                  run_name = "mixed_rates", conf = NULL),
                 "gene rate groups")

  # every sample gets the final-stage cumulative rate
  cesa = clear_gene_rates(cesa)
  cesa = set_gene_rates(cesa, rates = smp[! is.na(rate_late), .(gene, rate = rate_late)], missing_genes_take_nearest = TRUE)

  cesa = ces_variant_step(cesa, variants = kras_variants, stage_mut_prop = smp,
                          sample_index = sample_index, run_name = "step_test",
                          return_fit = TRUE, cores = 1, conf = 0.95)

  step_out = cesa@selection_results$step_test
  expect_true(all(c("si_early", "si_late", "loglikelihood") %in% names(step_out)))
  expect_true(all(c("ci_low_95_si_early", "ci_high_95_si_early",
                    "ci_low_95_si_late", "ci_high_95_si_late") %in% names(step_out)))
  expect_false("included_with_variant" %in% names(step_out))
  expect_true(all(step_out$si_early >= 0.001 & step_out$si_late >= 0.001))
  expect_true(!is.null(attr(step_out, "fit")))
  expect_true(all(c("selection_intensity_constant", "loglikelihood_constant") %in% names(step_out)))
  # constrained model is nested, so it can't fit better (beyond optimizer tolerance)
  expect_true(all(step_out$loglikelihood_constant <= step_out$loglikelihood + 1e-6))

  # requires exactly one gene per variant (or an unambiguous, valid stage_mut_prop row for it)
  expect_error(
    ces_variant_step(cesa, variants = select_variants(cesa, min_freq = 2), stage_mut_prop = smp,
                     sample_index = sample_index, run_name = "should_fail", cores = 1),
    "exactly.*one gene|not present in stage_mut_prop|NA stage_mut_prop values"
  )

  # LRT against the built-in constrained fit (default)
  lrt = step_selection_LRT(cesa, step_run_name = "step_test")
  expect_equal(lrt$loglik_simple, step_out$loglikelihood_constant)

  # The constrained fit is the same model as ces_variant() when each sample carries its own
  # stage's cumulative rate; check that, and the simple_run_name path of step_selection_LRT().
  nested = clear_gene_rates(cesa)
  nested = set_gene_rates(nested, rates = smp[! is.na(rate_early), .(gene, rate = rate_early)], samples = grpA, missing_genes_take_nearest = TRUE)
  nested = set_gene_rates(nested, rates = smp[! is.na(rate_late), .(gene, rate = rate_late)], samples = grpB, missing_genes_take_nearest = TRUE)
  nested = ces_variant(nested, variants = kras_variants, run_name = "simple_test", return_fit = TRUE, cores = 1)
  simple_out = nested@selection_results$simple_test
  expect_equal(simple_out$loglikelihood, step_out$loglikelihood_constant, tolerance = 1e-4)
  expect_equal(simple_out$selection_intensity, step_out$selection_intensity_constant, tolerance = 1e-3)
  lrt_ext = step_selection_LRT(nested, step_run_name = "step_test", simple_run_name = "simple_test")
  expect_equal(lrt_ext$p_value, lrt$p_value, tolerance = 1e-3)
  cesa = nested
  expect_setequal(names(lrt), c("variant_name", "loglik_step", "loglik_simple", "df", "LRT_stat", "p_value"))
  expect_true(all(lrt$df == 1))
  expect_true(all(lrt$p_value >= 0 & lrt$p_value <= 1))

  expect_error(step_selection_LRT(cesa, step_run_name = "not_a_run", simple_run_name = "simple_test"),
              "not found")

  # plot_effects_step(), skipped if ggplot2 unavailable
  skip_if_not_installed("ggplot2")
  p = plot_effects_step(step_out, group_by = "variant_name", lrt = lrt)
  expect_s3_class(p, "ggplot")
})
