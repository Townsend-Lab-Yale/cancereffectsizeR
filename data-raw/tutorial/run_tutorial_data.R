# load_all("~/cancereffectsizeR")
library(ggplot2)
library(ggrepel)


## Keep this script updated to match the quickstart

tcga_maf_file = 'TCGA.BRCA.mutect.995c0111-d90b-4140-bee7-3845436c3b42.DR-10.0.somatic.maf.gz'
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
                covered_regions = tgs_coverage, covered_regions_name = 'top_genes', covered_regions_padding = 1)

# Infer trinculeotide-context-specific relative rates of SNV mutation from
# a mutational signature analysis (leaving out signatures not found in BRCA)
signature_exclusions = suggest_cosmic_signatures_to_remove(cancer_type = 'BRCA', treatment_naive = TRUE)
cesa = trinuc_mutation_rates(cesa = cesa, signature_set = ces.refset.hg38$signatures$COSMIC_v3.2, 
                             signature_exclusions = signature_exclusions)

# Save snv_counts for visualization example
snv_count_file = paste0(system.file("tutorial/", package = 'cancereffectsizeR'), '/BRCA_snv_counts.rds')
sample_file = paste0(system.file("tutorial/", package = 'cancereffectsizeR'), '/BRCA_cesa_samples.rds')
saveRDS(cesa$mutational_signatures$snv_counts, snv_count_file)
saveRDS(cesa$samples, sample_file)

# Estimate neutral gene mutation rates using dNdScv, with tissue-specific mutation rate covariates.
cesa = gene_mutation_rates(cesa, covariates = ces.refset.hg38$covariates$breast)

# Infer scaled selection coefficients under the default model of clonal selection.
# By default, inference is restricted to recurrent mutations.
cesa = ces_variant(cesa)

# Copy output, merge in variant annotations, and view top variants
selection = cesa$selection[[1]]
selection = selection[cesa$variants, on = 'variant_id', nomatch = NULL]
top = selection[order(-selection_intensity)][1:10]

# Quickstart plot of top effects (also used for app note figure, panel B)
top = top[order(selection_intensity)] # will plot lowest to highest (left to right)

top_output = paste0(system.file("inst/tutorial/", package = 'cancereffectsizeR'), '/top_BRCA_effects.txt')
fwrite(top, top_output, sep = "\t")

top[, display_levels := factor(display_name, levels = display_name, ordered = T)]
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


# SNV counts
## Get pre-processed CESAnalysis 
met_tgs_maf = fread("")

