# Generates the pre-computed step-specific selection output shown in vignettes/step_specific_selection.Rmd
# and in the "Step-specific selection" section of vignettes/cancereffectsizeR.Rmd, saved to
# inst/step_selection_tutorial/.
#
# Data: esophageal squamous cell carcinoma (ESCC) analysis, normal esophageal epithelium ("Pre")
# vs. primary ESCC ("Pri"), with gene-level compound variants. The input CESAnalysis and
# CompoundVariantSet were built by the ESCC analysis (ESCC_step_epistasis/analysis/run_analysis.R)
# and are not distributed with the package. Only results for NOTCH1, TP53, and PIK3CA are saved.
# Needs ~100 GB of memory.
suppressMessages(library(ces.refset.hg19))
devtools::load_all("/gpfs/gibbs/project/townsend/kag227/escc_notch/cancereffectsizeR")
library(ggplot2)

in_dir = "/gpfs/gibbs/project/townsend/kag227/escc_notch/.step_selection_slurm/"
out_dir = "/gpfs/gibbs/project/townsend/kag227/escc_notch/cancereffectsizeR/inst/step_selection_tutorial"
genes = c("NOTCH1", "TP53", "PIK3CA")

cesa = readRDS(paste0(in_dir, "cesa_ready.rds"))
if (is.null(.ces_ref_data[[cesa@ref_key]])) .ces_ref_data[[cesa@ref_key]] = preload_ref_data(system.file("refset", package = cesa@ref_key))
compound = readRDS(paste0(in_dir, "compound.rds"))
compound = compound[which(compound@compounds$compound_name %in% genes)]

sample_index = assign_stage_index(cesa, stage_col = "Pre_or_Pri", stage_order = c("Pre", "Pri"))

cesa = clear_gene_rates(cesa)
cesa = gene_mutation_rates(cesa, covariates = "ESCA", save_all_dndscv_output = TRUE,
                           samples = cesa$samples[Pre_or_Pri == "Pre"])
cesa = gene_mutation_rates(cesa, covariates = "ESCA", save_all_dndscv_output = TRUE,
                           samples = cesa$samples[Pre_or_Pri == "Pri"])
smp = stage_mutation_proportions(cesa, rate_cols = c("rate_grp_1", "rate_grp_2"),
                                 stage_names = c("Pre", "Pri"))

cesa = clear_gene_rates(cesa)
cesa = set_gene_rates(cesa, rates = smp[, .(gene, rate = rate_Pri)], missing_genes_take_nearest = TRUE)
cesa = ces_variant_step(cesa, variants = compound, stage_mut_prop = smp, sample_index = sample_index,
                        run_name = "step_effects", conf = 0.95)
step_effects = cesa@selection_results$step_effects
lrt = step_selection_LRT(cesa, step_run_name = "step_effects")

saveRDS(smp[gene %in% genes], file.path(out_dir, "escc_stage_mut_prop.rds"))
saveRDS(step_effects, file.path(out_dir, "escc_step_effects.rds"))
saveRDS(lrt, file.path(out_dir, "escc_step_lrt.rds"))

p = plot_effects_step(step_effects, group_by = "variant_name", lrt = lrt, stage_order = c("Pre", "Pri"),
                      stage_labels = c(Pre = "Normal", Pri = "Tumor"))
ggsave(file.path(out_dir, "escc_step_plot.png"), p, width = 8, height = 4, dpi = 150)

print(smp[gene %in% genes])
print(step_effects)
print(lrt)
