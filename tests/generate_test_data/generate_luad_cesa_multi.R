prev_dir = setwd(system.file("tests/test_data/", package = "cancereffectsizeR"))
luad = CESAnalysis(genome = "hg19", progression_order = 1:4)
luad = load_maf(luad, maf = "luad.hg19.maf.txt", sample_col = "sample_id", 
                tumor_allele_col = "Tumor_Seq_Allele2", progression_col = "fake_stage")
luad = calc_baseline_mutation_rates(luad, covariate_file = "lung_pca")
saveRDS(luad, "cesa_for_snv_multi.rds")
test_genes = c("TTN", "KRAS", "RYR2", "EGFR", "TP53", "ASXL3","IFITM2")
luad = ces_snv(luad, genes = test_genes)
saveRDS(luad@selection_results, "multi_stage_snv_results.rds")


# repeat with subset of data for dndscv testing 
dndscv_samples = c("sample-1", "sample-106", "sample-108", "sample-11", "sample-31", "sample-33", "sample-35", 
                   "sample-40", "sample-46", "sample-6", "sample-67", "sample-68", "sample-7", "sample-71", 
                   "sample-73", "sample-74", "sample-77", "sample-82", "sample-83", "sample-95", "sample-99")
maf_for_dndscv = data.table::fread("luad.hg19.maf.txt")
maf_for_dndscv = maf_for_dndscv[sample_id %in% dndscv_samples]
for_dndscv = load_maf(cesa = CESAnalysis(genome="hg19", progression_order = 1:4), maf = maf_for_dndscv, sample_col = "sample_id",
                      tumor_allele_col = "Tumor_Seq_Allele2", progression_col = "fake_stage")
for_dndscv = trinucleotide_mutation_weights(for_dndscv)
saveRDS(for_dndscv, "cesa_for_multi_dndscv.rds")

# long tests will actually run dNdScv; short tests will just make sure internal preprocess/postprocess functions behave as expected
dndscv_input = cancereffectsizeR:::dndscv_preprocess(cesa = for_dndscv, covariate_file = "lung_pca")
saveRDS(dndscv_input, "dndscv_input_multi.rds")
dndscv_raw_output = lapply(dndscv_input, function(x) do.call(dndscv::dndscv, x))

# a few attributes are huge (>1 GB); drop these
dndscv_raw_output = lapply(dndscv_raw_output, function(x) { x$nbreg$terms = NULL; x$nbreg$model = NULL; x$poissmodel = NULL; return(x)})
saveRDS(dndscv_raw_output, "dndscv_raw_output_multi.rds")
dndscv_out = dndscv_postprocess(cesa = for_dndscv, dndscv_raw_output = dndscv_raw_output)


sel_cv = lapply(dndscv_out@dndscv_out_list, function(x) x$sel_cv)
saveRDS(sel_cv, "sel_cv_multi.rds")
saveRDS(dndscv_out@mutrates_list, "mutrates_multi.rds")
saveRDS(dndscv_out, "multi-stage-dndscv_pre-anno.rds") # for quick testing of annotation function 
anno_out = annotate_gene_maf(dndscv_out)
saveRDS(anno_out@annotated.snv.maf, "multi_annotated_maf_df.rds")

setwd(prev_dir)
