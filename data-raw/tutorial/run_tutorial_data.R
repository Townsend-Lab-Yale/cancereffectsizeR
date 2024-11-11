library(ggplot2)
library(devtools)
setwd("~/cancereffectsizeR/")
load_all()

tutorial_dir = system.file("tutorial", package = 'cancereffectsizeR') # verify that it's dev directory, not installed package
setwd(tutorial_dir)
# quickstart run (uses LUAD instead of BRCA)
tcga_maf_file = 'TCGA-LUAD.maf.gz'
if (! file.exists(tcga_maf_file)) {
  get_TCGA_project_MAF(project = 'LUAD', filename = tcga_maf_file)
}

# Prepare data
maf = preload_maf(maf = tcga_maf_file, refset = "ces.refset.hg38")

# Create cancereffectsizeR analysis and load data
cesa = CESAnalysis(refset = "ces.refset.hg38")
cesa = load_maf(cesa = cesa, maf = maf)

# Infer trinculeotide-context-specific relative rates of SNV mutation from
# a mutational signature analysis (leaving out signatures not found in LUAD)
signature_exclusions = suggest_cosmic_signature_exclusions(cancer_type = 'LUAD', treatment_naive = TRUE)
cesa = trinuc_mutation_rates(cesa = cesa, signature_set = ces.refset.hg38$signatures$COSMIC_v3.4,
                             signature_exclusions = signature_exclusions)

# Estimate neutral gene mutation rates using dNdScv, with tissue-specific mutation rate covariates.
cesa = gene_mutation_rates(cesa, covariates = ces.refset.hg38$covariates$lung)

# Infer scaled selection coefficients under the default model of clonal selection.
# By default, inference is restricted to recurrent mutations.
cesa = ces_variant(cesa, run_name = 'example')

# Visualize top-effect variants
p = plot_effects(effects = cesa$selection$example, color_by = '#DB382D', topn = 20)

# Attribute effects to mutational signatures
mut_effects = mutational_signature_effects(cesa, cesa$selection$example)

# Plot a comparison of how signatures contribute to mutation vs. selection
sig_effect_plot = plot_signature_effects(mut_effects, viridis_option = 'F', num_sig_groups = 5)

# Above, p = ... instead of just ggplot(...)
luad_plot_file = paste0(tutorial_dir, '/top_LUAD_effects.rds')
saveRDS(p, luad_plot_file)
luad_sig_effect_plot_file = paste0(tutorial_dir, '/LUAD_sig_effects.rds')
saveRDS(sig_effect_plot, luad_sig_effect_plot_file)
save_cesa(cesa, 'luad_quickstart.rds') # Note: Listed in .gitignore


# rest of tutorial uses BRCA
tcga_maf_file = 'TCGA-BRCA.maf.gz'
if (! file.exists(tcga_maf_file)) {
  get_TCGA_project_MAF(project = 'BRCA', filename = tcga_maf_file)
}


tcga_clinical = fread(system.file("tutorial/TCGA_BRCA_clinical.txt", package = "cancereffectsizeR"))
setnames(tcga_clinical, "patient_id", "Unique_Patient_Identifier") # change column name

tcga_maf = preload_maf(maf = tcga_maf_file, refset = "ces.refset.hg38")

# Prepare TGS data
tgs_maf_file = system.file('tutorial/metastatic_breast_2021_hg38.maf', package = 'cancereffectsizeR')
tgs_maf = preload_maf(maf = tgs_maf_file, refset = 'ces.refset.hg38')

# Create cancereffectsizeR analysis and load data
cesa = CESAnalysis(refset = "ces.refset.hg38")
cesa = load_maf(cesa = cesa, maf = tcga_maf, maf_name = 'BRCA')
cesa = load_sample_data(cesa, tcga_clinical)


top_tgs_genes = c("TP53", "PIK3CA", "ESR1","CDH1","GATA3","KMT2C",
                      "MAP3K1","AKT1","ARID1A","FOXA1","TBX3","PTEN")
tgs_coverage = ces.refset.hg38$gr_genes[ces.refset.hg38$gr_genes$gene %in% top_tgs_genes]

tgs_maf$pM = 'M1' # all samples metastatic
cesa = load_maf(cesa, maf = tgs_maf, sample_data_cols = 'pM', maf_name = 'MBC', coverage = 'targeted',
                covered_regions = tgs_coverage, covered_regions_name = 'top_genes', covered_regions_padding = 10)


# Infer trinculeotide-context-specific relative rates of SNV mutation from
# a mutational signature analysis (leaving out signatures not found in BRCA)
signature_exclusions = suggest_cosmic_signature_exclusions(cancer_type = 'BRCA', treatment_naive = TRUE)
cesa = trinuc_mutation_rates(cesa = cesa, signature_set = ces.refset.hg38$signatures$COSMIC_v3.4, 
                             signature_exclusions = signature_exclusions)

# Estimate neutral gene mutation rates using dNdScv, with tissue-specific mutation rate covariates.
cesa = gene_mutation_rates(cesa, covariates = ces.refset.hg38$covariates$breast)
saveRDS(cesa$gene_rates, paste0(tutorial_dir, '/BRCA_cesa_gene_rates.rds'))

dndscv_subset = cesa$dNdScv_results[[1]][qallsubs_cv < .05][order(qallsubs_cv)]
saveRDS(list(rate_grp_1 = dndscv_subset), paste0(tutorial_dir, '/BRCA_dndscv_out.rds'))

## Mutation rate example
variants_to_check = cesa$variants[order(-maf_prevalence), variant_id][1:3]
samples_to_check = c('TCGA-A2-A3Y0', 'TCGA-XX-A89A', 'P-0000224')
site_rates = baseline_mutation_rates(cesa = cesa, variant_ids = variants_to_check, samples = samples_to_check)
saveRDS(site_rates,  paste0(tutorial_dir, '/BRCA_site_rates_example.rds'))

# Infer scaled selection coefficients under the default model of clonal selection.
# By default, inference is restricted to recurrent mutations.
cesa = ces_variant(cesa, run_name = 'recurrents')

p = plot_effects(cesa$selection$recurrents)

brca_plot_file = paste0(tutorial_dir, '/top_BRCA_effects.rds')
saveRDS(p, brca_plot_file)

p2 = plot_effects(cesa$selection$recurrents, group_by = 'gene', topn = 10,
             label_individual_variants = FALSE)

brca_p2_file = paste0(tutorial_dir, '/BRCA_effects_in_top_genes.rds')
saveRDS(p2, brca_p2_file)

## Sequential 
# Take variants that appear at least twice in pM-annotated data
# We get IDs from the previously generated counts, and then take the corresponding entries
# from the CESAnalysis variants table.

# # (counts_by_M produced earlier in tutorial)
# counts_by_M = variant_counts(cesa = cesa, variant_ids = cesa$variants[maf_prevalence > 1, variant_id],
#                              by = 'pM')
# variants_for_sequential = counts_by_M[M0_prevalence + M1_prevalence > 2, variant_id]
# variants_for_sequential = cesa$variants[variants_for_sequential, on = 'variant_id']
# 
# 
# cesa = ces_variant(cesa, variants = variants_for_sequential, model = 'sequential', run_name = 'sequential', 
#                    ordering_col = 'pM', ordering = c('M0', 'M1'))
# 
# # Assess the same variants in the same samples using the single model
# cesa = ces_variant(cesa, variants = variants_for_sequential, model = 'basic', run_name = 'for_sequential_compare',
#                    samples = cesa$samples[!is.na(pM)])
# 
# combined_results = merge.data.table(cesa$selection$sequential, 
#                                     cesa$selection$for_sequential_compare, 
#                                     all.x = TRUE, all.y = FALSE, 
#                                     by = c('variant_id', 'variant_name', 'variant_type'), suffixes = c('.sequential', '.single'))
# 
# # Likelihood ratio test
# combined_results[, chisquared := -2 * (loglikelihood.single - loglikelihood.sequential)]
# combined_results[, p := pchisq(chisquared, df = 1, lower.tail = F)]
# 
# # Prep summary output for printing. Not shown here, but all of these variants are covered
# # in both of our data sources.
# for_print = combined_results[, .(variant_name, si_single = selection_intensity, si_M0, si_M1, p,
#                                  M0_count = included_with_variant_M0, M1_count = included_with_variant_M1)]
# sequential_signif_output = for_print[p < .05][order(p)]
# 
# web_table = copy(sequential_signif_output)
# web_table[, c("si_single", "si_M0", "si_M1", "p") := lapply(.SD, signif, 1), 
#           .SDcols = c("si_single", "si_M0", "si_M1", "p")]
# web_table[p>=.001, char_p := format(round(p, 3))]
# web_table[p < .001, char_p := format(signif(p, 1))]
# col_order = setdiff(names(web_table), "char_p")
# web_table[, p := NULL]
# setnames(web_table, "char_p", "p")
# setcolorder(web_table, col_order)
# sequential_output_file = paste0(tutorial_dir, '/sequential_signif_output.rds')
# saveRDS(web_table, sequential_output_file)

# Epistasis

## Variant epistasis
# Start by pulling full variant IDs (with protein identifier) from variants table
group1 = cesa$variants[c("PIK3CA E545K", "AKT1 E17K"), variant_id, on = 'variant_name']
group2 = cesa$variants[c("PIK3CA E545K", "PIK3CA E542K"), variant_id, on = 'variant_name']
cesa = ces_epistasis(cesa = cesa, variants = list(group1, group2), conf = .95, run_name = 'variant_epistasis_example')

variant_ep_output = paste0(tutorial_dir, '/variant_ep_example.rds')
saveRDS(cesa$epistasis$variant_epistasis_example, variant_ep_output)

## (Compound) variant epistasis
top_PIK3CA = cesa$variants[gene == 'PIK3CA' & maf_prevalence > 1]
top_akt1 = cesa$variants[variant_name == 'AKT1 E17K']
for_compound = rbind(top_PIK3CA, top_akt1)

# See define_compound_variants() documentation for details on arguments
comp = define_compound_variants(cesa = cesa, variant_table = for_compound, by = "gene", merge_distance = Inf)
cesa = ces_epistasis(cesa = cesa, variants = comp, run_name = "AKT1_E17K_vs_PIK3CA")
comp_ep_output = paste0(tutorial_dir, '/comp_variant_ep.rds')
saveRDS(cesa$epistasis$AKT1_E17K_vs_PIK3CA, comp_ep_output) 

## Gene epistasis
genes = c("AKT1", "PIK3CA", "TP53")

# Get consensus covered regions
combined_coverage = intersect(cesa$coverage_ranges$exome$`exome+`, cesa$coverage_ranges$targeted$top_genes)

# Get variants in the genes of interest that have sequencing coverage in all samples
variants = select_variants(cesa, genes = genes, gr = combined_coverage)

cesa = ces_gene_epistasis(cesa = cesa, genes = genes, variants = variants, run_name = "gene_epistasis_example")
gene_ep_output = paste0(tutorial_dir, '/gene_ep_example.rds')
saveRDS(cesa$epistasis$gene_epistasis_example, gene_ep_output)


results = cesa$epistasis$gene_epistasis_example
epistasis_plot = plot_epistasis(results)[[1]]
brca_epistasis_file = paste0(tutorial_dir, '/BRCA_epistasis_example.rds')
saveRDS(epistasis_plot, brca_epistasis_file)

## Save CESAnalysis for reference/revisions
save_cesa(cesa, 'brca_tutorial_cesa.rds')

## Ratio-style epistasis plot

# results = results[, .(v1 = variant_A, v2 = variant_B, ces_A0, ces_B0, ces_A_on_B, 
#                                 ces_B_on_A, p_A_change, p_B_change, p_epistasis)]
# 
# # By change, we mean fold-change of selection on mutant background over wildtype background
# results[, change_in_v2 := ces_B_on_A / ces_B0]
# results[, change_in_v1 := ces_A_on_B /ces_A0]
# 
# # Put in desired order for display
# results = results[, pairname := paste(v1, v2, sep = '.')]
# 
# results[, x := 1:.N]
# results[, v1_x := x - .2]
# results[, v2_x := x + .2]
# results[, alpha := .6]
# results[p_epistasis < .05, alpha := 1]
# 
# results[, v1_signif := ""]
# results[p_A_change < .05, v1_signif := "*"]
# results[p_A_change < .01, v1_signif := "**"]
# results[p_A_change < .001, v1_signif := "***"]
# results[, v1_signif_y := change_in_v1 + (.13 * sign(change_in_v1 - 1))]
# 
# results[, v2_signif := ""]
# results[p_B_change < .05, v2_signif := "*"]
# results[p_B_change < .01, v2_signif := "**"]
# results[p_B_change < .001, v2_signif := "***"]
# results[, v2_signif_y := change_in_v2 + (.13 * sign(change_in_v2 - 1))]
# 
# x_labels = unlist(S4Vectors::zipup(results$v1, results$v2))
# x_label_pos = unlist(S4Vectors::zipup(results$v1_x, results$v2_x))
# 
# # Have to get fancy to depict significance nicely in legend.
# draw_signif_key = function(data, params, size) {
#   grobTree(rectGrob(x = .25, y = .5, width = .5, height = 1, 
#                     gp = gpar(col = NA, fill = alpha('plum4', data$alpha), lty = data$linetype)),
#            rectGrob(x = .75, y = .5, width = .5, height = 1, 
#                     gp = gpar(col = NA, fill = alpha('sandybrown', data$alpha), lty = data$linetype)))
# }
# 
# epistasis_plot = ggplot(data = results) + 
#   # Put in a reference line depicting no change in selection
#   geom_hline(yintercept = 1, color = 'darkgrey') +
#   
#   geom_rect(aes(xmin = v1_x - .2, xmax = v1_x + .2, ymin = 1, ymax = change_in_v1, fill = 'v1', alpha = alpha),
#             show.legend = c(alpha = FALSE, fill = TRUE)) +
#   geom_rect(aes(xmin = 1, xmax = 1, ymin = 0, ymax = 0, alpha = alpha), 
#             show.legend = c(alpha = TRUE, fill = FALSE), key_glyph = draw_signif_key) +
#   geom_text(aes(x = v1_x, y = v1_signif_y, label = v1_signif), size = 7) +
#   scale_alpha_identity(breaks = c(1, .6), labels = c('Significant', 'Not significant'),
#                        guide = guide_legend(title = 'Pairwise epistasis', override.aes = list(fill = 'sandybrown', alpha = c(1, .6)),
#                                             order = 1)) +
#   
#   geom_rect(aes(xmin = v2_x - .2, xmax = v2_x + .2, ymin = 1, ymax = change_in_v2, fill = 'v2', alpha = alpha),
#             show.legend = c(alpha = FALSE, fill = TRUE)) +
#   geom_text(aes(x = v2_x, y = v2_signif_y, label = v2_signif), size = 7) + 
#   
#   # Build legend
#   scale_fill_manual(name = 'Ratio of selection', 
#                     breaks = c('v1', 'v2'), 
#                     labels = list(expression(frac('gene 1 on mutated gene 2', 'gene 1 on wildtype gene 2')), 
#                                   expression(frac('gene 2 on mutated gene 1', 'gene 2 on wildtype gene 1'))),
#                     values = c('v1' = 'plum4', 'v2' = 'sandybrown'),
#                     guide = guide_legend(label.theme = element_text(size = 6.5))) +
#   scale_x_continuous(breaks = x_label_pos, labels = x_labels) + 
#   scale_y_continuous(breaks = seq(from = 0, to = 3, by = .25)) +
#   xlab('Gene pair') + ylab('Ratio of selection coefficients') + 
#   theme_classic() +
#   theme(legend.position = 'bottom', legend.title = element_text(size = 10), 
#         axis.ticks.length.x = unit(0, 'cm'))
#


