# Generates the pre-computed example output shown in vignettes/step_specific_selection.Rmd,
# saved to inst/step_selection_tutorial/. Run from the package root.
#
# The example is synthetic: the package's hg38 test CESAnalysis, split arbitrarily into "Normal"
# and "Tumor" halves, with made-up gene mutation rates. It demonstrates workflow and output format
# only.
devtools::load_all(".")
library(ggplot2)

out_dir = "inst/step_selection_tutorial"
dir.create(out_dir, showWarnings = FALSE)

cesa = load_cesa("tests/test_data/cesa_hg38_for_test.rds")
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
  group_name = c(rep("Normal", length(grpA)), rep("Tumor", length(grpB)))
)
setkey(sample_index, "Unique_Patient_Identifier")

# nearest-gene imputation of the genes without made-up rates produces many non-monotonic genes
smp = suppressWarnings(stage_mutation_proportions(cesa, rate_cols = c("rate_grp_1", "rate_grp_2"),
                                                  stage_names = c("Normal", "Tumor")))

kras_variants = select_variants(cesa, genes = "KRAS", min_freq = 2)
cesa = ces_variant_step(cesa, variants = kras_variants, stage_mut_prop = smp,
                        sample_index = sample_index, run_name = "step_effects", conf = 0.95)
step_effects = cesa@selection_results$step_effects

cesa = ces_variant(cesa, variants = kras_variants, run_name = "simple_effects", return_fit = TRUE)
lrt = step_selection_LRT(cesa, step_run_name = "step_effects", simple_run_name = "simple_effects")

saveRDS(smp[gene %in% genes], file.path(out_dir, "kras_stage_mut_prop.rds"))
saveRDS(step_effects, file.path(out_dir, "kras_step_effects.rds"))
saveRDS(lrt, file.path(out_dir, "kras_step_lrt.rds"))

p = plot_effects_step(step_effects, group_by = "variant_name", lrt = lrt, stage_order = c("Normal", "Tumor"))
ggsave(file.path(out_dir, "kras_step_plot.png"), p, width = 8, height = 5, dpi = 150)
