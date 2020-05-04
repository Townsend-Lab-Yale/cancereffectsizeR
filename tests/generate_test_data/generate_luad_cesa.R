prev_dir = setwd(system.file("tests/test_data/", package = "cancereffectsizeR"))

# read in the MAF used for all the testing
luad = load_maf(cesa = CESAnalysis(genome="hg19"), maf = "luad.hg19.maf.txt", sample_col = "sample_id")

# use the trinucleotide weight data saved in generate_trinuc_data (seldom a need to change that data)
weights = readRDS("luad_trinucleotide_mutation_weights.rds")
luad = set_trinuc_rates(luad, weights$trinuc_proportion_matrix)
luad = gene_mutation_rates(luad, covariate_file = "lung_pca")
luad = annotate_variants(luad)

saveRDS(luad, "cesa_for_snv.rds")

# save results to serve as expected test output
test_genes = c("TTN", "EGFR", "ASXL3", "KRAS", "RYR2", "USH2A", "CSMD3", "TP53", "CSMD1", "LRP1B", 
               "ZFHX4", "FAT3", "CNTNAP5", "PCDH15", "NEB", "RYR3", "DMD", "KATNAL1", 
               "OR13H1", "KSR1")
luad = ces_snv(luad, genes = test_genes)
saveRDS(luad@selection_results, "single_stage_snv_results.rds")
luad = ces_gene_epistasis(luad, genes = c("EGFR", "KRAS", "TP53"), return_all_opm_output = T)
saveRDS(luad@gene_epistasis_results, "epistasis_results.rds")
saveRDS(luad@advanced$opm_output, "epistasis_opm.rds")


# repeat with smaller data set for dndscv testing
for_dndscv = load_maf(cesa = CESAnalysis(genome="hg19"), maf = "luad.hg19.maf.txt", sample_col = "sample_id")
dndscv_samples = c("sample-1", "sample-106", "sample-108", "sample-11", "sample-31", "sample-33", "sample-35", 
                   "sample-40", "sample-46", "sample-6", "sample-67", "sample-68", "sample-7", "sample-71", 
                   "sample-73", "sample-74", "sample-77", "sample-82", "sample-83", "sample-95", "sample-99")

maf_for_dndscv = luad@maf[Unique_Patient_Identifier %in% dndscv_samples]

for_dndscv = load_maf(cesa = CESAnalysis(genome="hg19"), maf = maf_for_dndscv, sample_col = "Unique_Patient_Identifier", tumor_allele_col = "Tumor_Allele")
for_dndscv = set_trinuc_rates(for_dndscv, weights$trinuc_proportion_matrix[dndscv_samples,])
saveRDS(for_dndscv, "cesa_for_dndscv_and_anno.rds")

# long tests will actually run dNdScv; short tests will just make sure internal preprocess/postprocess functions behave as expected
dndscv_input = cancereffectsizeR:::dndscv_preprocess(cesa = for_dndscv, covariate_file = "lung_pca")
saveRDS(dndscv_input, "dndscv_input_single.rds")
dndscv_raw_output = lapply(dndscv_input, function(x) do.call(dndscv::dndscv, x))

# a few attributes are huge (>1 GB); drop these
dndscv_raw_output = lapply(dndscv_raw_output, function(x) { x$nbreg$terms = NULL; x$nbreg$model = NULL; x$poissmodel = NULL; return(x)})
saveRDS(dndscv_raw_output, "dndscv_raw_output_single.rds")
dndscv_out = dndscv_postprocess(cesa = for_dndscv, dndscv_raw_output = dndscv_raw_output)

saveRDS(dndscv_out@dndscv_out_list[[1]]$sel_cv, "sel_cv.rds")
saveRDS(dndscv_out@mutrates_list[[1]], "mutrates.rds")
saveRDS(dndscv_out, "single-dndscv_pre-anno.rds")
anno_out = annotate_variants(dndscv_out)
saveRDS(anno_out@maf, "annotated_maf_df.rds")

setwd(prev_dir)




