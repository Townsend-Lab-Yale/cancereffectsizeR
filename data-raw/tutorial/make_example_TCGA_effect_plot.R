# load_all("~/cancereffectsizeR")
library(ggplot2)
library(ggrepel)


## Keep this updated to essentially match the quickstart
tcga_maf_file = 'TCGA.BRCA.mutect.995c0111-d90b-4140-bee7-3845436c3b42.DR-10.0.somatic.maf.gz'
if (! file.exists(tcga_maf_file)) {
  download.file('https://api.gdc.cancer.gov/data/995c0111-d90b-4140-bee7-3845436c3b42', 
                destfile = tcga_maf_file)
}

maf = preload_maf(maf = tcga_maf_file, refset = "ces.refset.hg38")
maf = maf[germline_variant_site == F][repetitive_region == F | cosmic_site_tier %in% 1:3]

# Create cancereffectsizeR analysis and load data
cesa = CESAnalysis(refset = "ces.refset.hg38")
cesa = load_maf(cesa = cesa, maf = maf)

# Infer trinculeotide-context-specific relative rates of SNV mutation from
# a mutational signature analysis (leaving out signatures not found in BRCA)
signature_exclusions = suggest_cosmic_signatures_to_remove(cancer_type = 'BRCA', treatment_naive = TRUE)
cesa = trinuc_mutation_rates(cesa = cesa, signature_set = ces.refset.hg38$signatures$COSMIC_v3.2, 
                             signature_exclusions = signature_exclusions)

# Estimate neutral gene mutation rates using dNdScv, with tissue-specific mutation rate covariates.
cesa = gene_mutation_rates(cesa, covariates = ces.refset.hg38$covariates$breast)

# Infer scaled selection coefficients under the default model of clonal selection.
# By default, inference is restricted to recurrent mutations.
cesa = ces_variant(cesa)

# Copy output, merge in variant annotations, and view top variants
selection = cesa$selection[[1]]
selection = selection[cesa$variants, on = 'variant_id', nomatch = NULL]
top = selection[order(-selection_intensity)][1:10]

# Plot top effects (also used for app note figure, panel B)
top = top[order(selection_intensity)] # will plot lowest to highest (left to right)
top[, display_levels := factor(display_name, levels = display_name, ordered = T)]

top_output = paste0(system.file("inst/tutorial/", package = 'cancereffectsizeR'), '/top_BRCA_effects.txt')
fwrite(top, top_output, sep = "\t")

ggplot(top, aes(display_levels, selection_intensity)) + 
  geom_errorbar(aes(ymin = ci_low_95, ymax = ci_high_95), width = .2, color = 'darkgrey') +
  scale_x_discrete() + scale_y_log10() +
  geom_point(aes(color = as.factor(maf_prevalence)), size = 2.5) + 
  xlab("variant") + ylab(expression('cancer effect'~scriptstyle(~~(log[10])))) +
  ggtitle('Top cancer effects in TCGA BRCA') +
  guides(color = guide_legend(title = 'variant prevalence')) + 
  theme_minimal() + theme(legend.position = 'bottom', text = element_text(family = "Arial"),
                          axis.text.x = element_text(size = 8, family = "Verdana"))





## Alternate plot
# top[, display_levels := factor(display_name, levels = display_name, ordered = T)]
# ggplot(top, aes(display_levels, selection_intensity)) + 
#   geom_errorbar(aes(ymin = ci_low_95, ymax = ci_high_95), width = .2, color = 'darkgrey') +
#   scale_x_discrete() + scale_y_log10() +
#   geom_point(aes(color = as.factor(maf_prevalence)), size = 2.5) + 
#   xlab("variant") + ylab(expression('cancer effect'~scriptstyle(~~(log[10])))) +
#   ggtitle('Top cancer effects in TCGA BRCA') +
#   guides(color = guide_legend(title = 'variant prevalence')) + 
#   theme_minimal() + theme(legend.position = 'bottom', text = element_text(family = "Helvetica"),
#                           axis.text.x = element_text(size = 8, family = "Verdana"))

