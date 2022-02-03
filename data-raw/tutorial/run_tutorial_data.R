library(ggplot2)
library(devtools)
setwd("~/cancereffectsizeR/")
load_all()

tutorial_dir = system.file("tutorial", package = 'cancereffectsizeR') # verify that it's dev directory, not installed package

# quickstart run (uses LUAD instead of BRCA)
tcga_maf_file = 'inst/tutorial/TCGA.LUAD.mutect.0458c57f-316c-4a7c-9294-ccd11c97c2f9.DR-10.0.somatic.maf.gz'
if (! file.exists(tcga_maf_file)) {
  download.file('https://api.gdc.cancer.gov/data/0458c57f-316c-4a7c-9294-ccd11c97c2f9', 
                destfile = tcga_maf_file)
}

maf = preload_maf(maf = tcga_maf_file, refset = "ces.refset.hg38")
maf = maf[germline_variant_site == F][repetitive_region == F | cosmic_site_tier %in% 1:3]

# Create cancereffectsizeR analysis and load data
cesa = CESAnalysis(refset = "ces.refset.hg38")
cesa = load_maf(cesa = cesa, maf = maf)

# Infer trinculeotide-context-specific relative rates of SNV mutation from
# a mutational signature analysis (leaving out signatures not found in LUAD)
signature_exclusions = suggest_cosmic_signature_exclusions(cancer_type = 'LUAD', treatment_naive = TRUE)
cesa = trinuc_mutation_rates(cesa = cesa, signature_set = ces.refset.hg38$signatures$COSMIC_v3.2,
                             signature_exclusions = signature_exclusions)

# Estimate neutral gene mutation rates using dNdScv, with tissue-specific mutation rate covariates.
cesa = gene_mutation_rates(cesa, covariates = ces.refset.hg38$covariates$lung)

# Infer scaled selection coefficients under the default model of clonal selection.
# By default, inference is restricted to recurrent mutations.
cesa = ces_variant(cesa)

# Take top 15 variants, then sort lowest to highest (to plot left to right)
selection = cesa$selection[[1]]
top = selection[order(-selection_intensity)][1:15]
top = top[order(selection_intensity)]

# Plot top effects 
top[, display_name := gsub('_', "\n", variant_name)]
top[, display_levels := factor(display_name, levels = display_name, ordered = T)]
plot_title = 'Top cancer effects in TCGA LUAD'
breaks = unique(as.numeric(round(quantile(top$included_with_variant))))

p = ggplot(top, aes(x = display_levels, y = selection_intensity)) + 
  geom_errorbar(aes(ymin = ci_low_95, ymax = ci_high_95), width = .2, color = 'darkgrey') +
  geom_point(aes(color = included_with_variant), size = 3) + 
  scale_x_discrete() + scale_y_log10() + 
  scale_color_viridis_c(name = 'variant prevalence', guide = 'colorbar', trans = 'log10', 
                        option = 'plasma', breaks = breaks) +
  xlab(element_blank()) +
  ylab(expression('cancer effect'~scriptstyle(~~(log[10])))) +
  ggtitle(plot_title) +
  guides(color = guide_colourbar(ticks = FALSE)) + 
  theme_minimal() + 
  theme(text = element_text(family = "Verdana"),
        axis.title.x = element_text(size = 14), 
        axis.text.x = element_text(size = 8),
        legend.position = 'bottom',
        legend.title = element_text(size = 10),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

luad_plot_file = paste0(tutorial_dir, '/top_LUAD_effects.rds')
saveRDS(p, luad_plot_file)


# rest of tutorial uses BRCA
tcga_maf_file = 'inst/tutorial/TCGA.BRCA.mutect.995c0111-d90b-4140-bee7-3845436c3b42.DR-10.0.somatic.maf.gz'
if (! file.exists(tcga_maf_file)) {
  download.file('https://api.gdc.cancer.gov/data/995c0111-d90b-4140-bee7-3845436c3b42', 
                destfile = tcga_maf_file)
}

tcga_maf = preload_maf(maf = tcga_maf_file, refset = "ces.refset.hg38")
tcga_maf = tcga_maf[germline_variant_site == F][repetitive_region == F | cosmic_site_tier %in% 1:3]

tcga_clinical = fread(system.file("tutorial/TCGA_BRCA_clinical.txt", package = "cancereffectsizeR"))
tcga_maf[, patient_id := substr(Unique_Patient_Identifier, 1, 12)]

uniqueN(tcga_maf[, .(patient_id, Unique_Patient_Identifier)]) == uniqueN(tcga_maf$Unique_Patient_Identifier)

# We can now use patient_id as the unique identifier.
tcga_maf[, Unique_Patient_Identifier := patient_id]
tcga_maf[, patient_id := NULL] # remove redundant column
setnames(tcga_clinical, "patient_id", "Unique_Patient_Identifier") # change column name

# Create cancereffectsizeR analysis and load data
cesa = CESAnalysis(refset = "ces.refset.hg38")
cesa = load_maf(cesa = cesa, maf = tcga_maf)
cesa = load_sample_data(cesa, tcga_clinical)

# Load in TGS data
metastatic_tgs_maf_file = system.file("tutorial/metastatic_breast_2021_hg38.maf", package = "cancereffectsizeR")
metastatic_tgs_maf = fread(metastatic_tgs_maf_file)
metastatic_tgs_maf$pM = 'M1'

top_tgs_genes = c("TP53", "PIK3CA", "ESR1","CDH1","GATA3","KMT2C",
                      "MAP3K1","AKT1","ARID1A","FOXA1","TBX3","PTEN")
tgs_coverage = ces.refset.hg38$gr_genes[ces.refset.hg38$gr_genes$gene %in% top_tgs_genes]
cesa = load_maf(cesa, maf = metastatic_tgs_maf, sample_data_cols = 'pM', coverage = 'targeted',
                covered_regions = tgs_coverage, covered_regions_name = 'top_genes')

# Infer trinculeotide-context-specific relative rates of SNV mutation from
# a mutational signature analysis (leaving out signatures not found in BRCA)
signature_exclusions = suggest_cosmic_signature_exclusions(cancer_type = 'BRCA', treatment_naive = TRUE)
cesa = trinuc_mutation_rates(cesa = cesa, signature_set = ces.refset.hg38$signatures$COSMIC_v3.2, 
                             signature_exclusions = signature_exclusions)

# Save snv_counts for visualization example
snv_count_file = paste0(tutorial_dir, '/BRCA_snv_counts.rds')
sample_file = paste0(tutorial_dir, '/BRCA_cesa_samples.rds')
saveRDS(cesa$mutational_signatures$snv_counts, snv_count_file)
saveRDS(cesa$samples, sample_file)

# Estimate neutral gene mutation rates using dNdScv, with tissue-specific mutation rate covariates.
cesa = gene_mutation_rates(cesa, covariates = ces.refset.hg38$covariates$breast)
saveRDS(cesa$gene_rates, paste0(tutorial_dir, '/BRCA_cesa_gene_rates.rds'))

dndscv_subset = cesa$dNdScv_results[[1]][qallsubs_cv < .05, .SD[which.min(qallsubs_cv)], by = 'gene'][order(qallsubs_cv)]
saveRDS(list(rate_grp_1 = dndscv_subset), paste0(tutorial_dir, '/BRCA_dndscv_out.rds'))

## Mutation rate example
variants_to_check = cesa$variants[order(-maf_prevalence), variant_id][1:3]
samples_to_check = c('TCGA-A2-A3Y0', 'TCGA-XX-A89A', 'P-0000224')
site_rates = baseline_mutation_rates(cesa = cesa, variant_ids = variants_to_check, samples = samples_to_check)
saveRDS(site_rates,  paste0(tutorial_dir, '/BRCA_site_rates_example.rds'))


# Infer scaled selection coefficients under the default model of clonal selection.
# By default, inference is restricted to recurrent mutations.
cesa = ces_variant(cesa, run_name = 'recurrents')

# Extract selection results from CESAnalysis and take top variants for visualization
top = cesa$selection$recurrents
top = top[order(-selection_intensity)][1:20] # take top 20 by SI
top = top[order(selection_intensity)] # will plot lowest to highest (left to right)

# Make variant names pretty for use in plot labels
top[, display_name := gsub('_', "\n", variant_name)]
top[, display_levels := factor(display_name, levels = display_name, ordered = T)]

plot_title = 'Top cancer effects in breast carcinoma (CES tutorial data)'
breaks = unique(as.numeric(round(quantile(top$included_with_variant))))
n.dodge = 2 # can reduce to 1 if labels happen to still fit (e.g., if plotting fewer variants)
p = ggplot(top, aes(x = display_levels, y = selection_intensity)) + 
  geom_errorbar(aes(ymin = ci_low_95, ymax = ci_high_95), width = .2, color = 'darkgrey') +
  geom_point(aes(color = included_with_variant), size = 3) + 
  scale_x_discrete(guide = guide_axis(n.dodge = n.dodge)) + scale_y_log10() + 
  scale_color_viridis_c(name = 'variant prevalence', guide = 'colorbar', trans = 'log10', 
                        option = 'plasma', breaks = breaks) +
  xlab(element_blank()) +
  ylab(expression('cancer effect'~scriptstyle(~~(log[10])))) +
  ggtitle(plot_title) +
  guides(color = guide_colourbar(ticks = FALSE)) + 
  theme_minimal() + 
  theme(text = element_text(family = "Verdana"),
        axis.title.x = element_text(size = 14), 
        axis.text.x = element_text(size = 8),
        legend.position = 'bottom',
        legend.title = element_text(size = 10),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

brca_plot_file = paste0(tutorial_dir, '/top_BRCA_effects.rds')
saveRDS(p, brca_plot_file)

## Sequential 
# Take variants that appear at least twice in pM-annotated data
# We get IDs from the previously generated counts, and then take the corresponding entries
# from the CESAnalysis variants table.

# (counts_by_M produced earlier in tutorial)
counts_by_M = variant_counts(cesa = cesa, variant_ids = cesa$variants[maf_prevalence > 1, variant_id],
                             by = 'pM')
variants_for_sequential = counts_by_M[M0_prevalence + M1_prevalence > 2, variant_id]
variants_for_sequential = cesa$variants[variants_for_sequential, on = 'variant_id']


cesa = ces_variant(cesa, variants = variants_for_sequential, model = 'sequential', run_name = 'sequential', 
                   ordering_col = 'pM', ordering = c('M0', 'M1'))

# Assess the same variants in the same samples using the single model
cesa = ces_variant(cesa, variants = variants_for_sequential, model = 'basic', run_name = 'for_sequential_compare',
                   samples = cesa$samples[!is.na(pM)])

combined_results = merge.data.table(cesa$selection$sequential, 
                                    cesa$selection$for_sequential_compare, 
                                    all.x = TRUE, all.y = FALSE, 
                                    by = c('variant_id', 'variant_name', 'variant_type'), suffixes = c('.sequential', '.single'))

# Likelihood ratio test
combined_results[, chisquared := -2 * (loglikelihood.single - loglikelihood.sequential)]
combined_results[, p := pchisq(chisquared, df = 1, lower.tail = F)]

# Prep summary output for printing. Not shown here, but all of these variants are covered
# in both of our data sources.
for_print = combined_results[, .(variant_name, si_single = selection_intensity, si_M0, si_M1, p,
                                 M0_count = included_with_variant_M0, M1_count = included_with_variant_M1)]
sequential_signif_output = for_print[p < .05][order(p)]

web_table = copy(sequential_signif_output)
web_table[, c("si_single", "si_M0", "si_M1", "p") := lapply(.SD, signif, 1), 
          .SDcols = c("si_single", "si_M0", "si_M1", "p")]
web_table[p>=.001, char_p := format(round(p, 3))]
web_table[p < .001, char_p := format(signif(p, 1))]
col_order = setdiff(names(web_table), "char_p")
web_table[, p := NULL]
setnames(web_table, "char_p", "p")
setcolorder(web_table, col_order)
sequential_output_file = paste0(tutorial_dir, '/sequential_signif_output.rds')
saveRDS(web_table, sequential_output_file)

# Epistasis

## Variant epistasis
# Start by pulling full variant IDs (with protein identifier) from variants table
group1 = cesa$variants[c("PIK3CA_E545K", "AKT1_E17K"), variant_id, on = 'variant_name']
group2 = cesa$variants[c("PIK3CA_E545K", "PIK3CA_E542K"), variant_id, on = 'variant_name']
cesa = ces_epistasis(cesa = cesa, variants = list(group1, group2), conf = .95, run_name = 'variant_epistasis_example')

variant_ep_output = paste0(tutorial_dir, '/variant_ep_example.rds')
saveRDS(cesa$epistasis$variant_epistasis_example, variant_ep_output)

## (Compound) variant epistasis
top_PIK3CA = cesa$variants[gene == 'PIK3CA' & maf_prevalence > 1]
top_akt1 = cesa$variants[variant_name == 'AKT1_E17K']
for_compound = rbind(top_PIK3CA, top_akt1)

# See define_compound_variants() documentation for details on arguments
comp = define_compound_variants(cesa = cesa, variant_table = for_compound, by = "gene", merge_distance = Inf)
cesa = ces_epistasis(cesa = cesa, variants = comp, run_name = "AKT1_E17K_vs_PIK3CA")
comp_ep_output = paste0(tutorial_dir, '/comp_variant_ep.rds')
saveRDS(cesa$epistasis$AKT1_E17K_vs_PIK3CA, comp_ep_output) 

## Gene epistasis
genes = c("AKT1", "PIK3CA", "TP53")
variants = cesa$variants[variant_type == 'aac' & gene %in% genes & sapply(covered_in, length) == 2]
variants = variants[aa_ref != aa_alt | essential_splice == T]
cesa = ces_gene_epistasis(cesa = cesa, genes = genes, variants = variants, run_name = "gene_epistasis_example")
gene_ep_output = paste0(tutorial_dir, '/gene_ep_example.rds')
saveRDS(cesa$epistasis$gene_epistasis_example, gene_ep_output)


## Save CESAnalysis for reference/revisions
#save_cesa(cesa, 'brca_cesa.rds')

