library(ggplot2)
library(devtools)
setwd("~/cancereffectsizeR/")
load_all()


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

# Copy output, merge in variant annotations, and view top variants
selection = cesa$selection[[1]]
selection = selection[cesa$variants, on = 'variant_id', nomatch = NULL]


# Take top 15 variants, then sort lowest to highest (to plot left to right)
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

luad_plot_file = paste0(system.file("tutorial/", package = 'cancereffectsizeR'), '/top_LUAD_effects.rds')
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
snv_count_file = paste0(system.file("tutorial/", package = 'cancereffectsizeR'), '/BRCA_snv_counts.rds')
sample_file = paste0(system.file("tutorial/", package = 'cancereffectsizeR'), '/BRCA_cesa_samples.rds')
saveRDS(cesa$mutational_signatures$snv_counts, snv_count_file)
saveRDS(cesa$samples, sample_file)

# Estimate neutral gene mutation rates using dNdScv, with tissue-specific mutation rate covariates.
cesa = gene_mutation_rates(cesa, covariates = ces.refset.hg38$covariates$breast)
saveRDS(cesa$gene_rates, 
        paste0(system.file("tutorial/", package = 'cancereffectsizeR'), '/BRCA_cesa_gene_rates.rds'))

dndscv_subset = cesa$dNdScv_results[[1]][qallsubs_cv < .05, .SD[which.min(qallsubs_cv)], by = 'gene'][order(qallsubs_cv)]
saveRDS(list(rate_grp_1 = dndscv_subset), paste0(system.file("tutorial/", package = 'cancereffectsizeR'), '/BRCA_dndscv_out.rds'))

## Mutation rate example
variants_to_check = cesa$variants[order(-maf_prevalence), variant_id][1:3]
samples_to_check = c('TCGA-A2-A3Y0', 'TCGA-XX-A89A', 'P-0000224')
site_rates = baseline_mutation_rates(cesa = cesa, variant_ids = variants_to_check, samples = samples_to_check)
saveRDS(site_rates,  paste0(system.file("tutorial/", package = 'cancereffectsizeR'), '/BRCA_site_rates_example.rds'))


# Infer scaled selection coefficients under the default model of clonal selection.
# By default, inference is restricted to recurrent mutations.
cesa = ces_variant(cesa, run_name = 'recurrents')

# Copy output, merge in variant annotations, and view top variants
top = cesa$selection$recurrents
top = top[cesa$variants, on = 'variant_id', nomatch = NULL]
top = top[order(-selection_intensity)][1:15]
top = top[order(selection_intensity)] # will plot lowest to highest (left to right)

# Plot top effects 
top[, display_name := gsub('_', "\n", variant_name)]
top[, display_levels := factor(display_name, levels = display_name, ordered = T)]


plot_title = 'Top cancer effects in breast carcinoma (CES tutorial data)'
breaks = as.numeric(round(quantile(top$included_with_variant, c(0, .33, .66, 1))))
ggplot(top, aes(x = display_levels, y = selection_intensity)) + 
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


variants_for_sequential = cesa$variants[maf_prevalence > 2]
cesa = ces_variant(cesa, variants = variants_for_sequential, model = 'sequential', run_name = 'sequential', 
                   ordering_col = 'pM', ordering = c('M0', 'M1'))

of_interest = cesa$selection$sequential[ci_high_95_si_M0 < ci_low_95_si_M1 | ci_high_95_si_M1 < ci_low_95_si_M0]

variants_for_compound = cesa$variants[gene %in% c('ESR1', 'PIK3CA', 'TP53')][sapply(covered_in, length) == 2]
variants_for_compound[aa_ref == aa_alt & essential_splice == F, table(maf_prevalence)]
variants_for_compound = variants_for_compound[aa_ref != aa_alt | essential_splice == T | variant_type == 'snv']
comp = define_compound_variants(cesa, variant_table = variants_for_compound, by = 'gene', merge_distance = Inf)

cesa = ces_variant(cesa, variants = comp, model = 'sequential', run_name = 'comp-sequential', ordering_col = 'pM', ordering = c('M0', 'M1'))


# Let's look at c("AKT1", 


cesa = ces_gene_epistasis(cesa, genes = c('ESR1', 'GATA3', 'FOXA1'))


top_output = paste0(system.file("inst/tutorial/", package = 'cancereffectsizeR'), '/top_TCGA_effects_for_quickstart.txt')
fwrite(top, top_output, sep = "\t")


