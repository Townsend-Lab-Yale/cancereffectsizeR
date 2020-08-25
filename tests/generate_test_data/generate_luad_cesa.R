prev_dir = setwd(system.file("tests/test_data/", package = "cancereffectsizeR"))

# read in the MAF used for all the testing
luad = load_maf(cesa = CESAnalysis(genome="hg19"), maf = "luad.hg19.maf.txt", sample_col = "sample_id")

# use the trinucleotide weight data saved in generate_trinuc_data (seldom a need to change that data)
trinuc_rates = fread("luad_hg19_trinuc_rates.txt")
luad = set_trinuc_rates(luad, trinuc_rates = trinuc_rates)


# Save gene_mutation_rates intermediary stuff so that tests don't need to run dNdScv
# This will not be the same as running gene_mutation_rates directly due to dNdScv gr_genes stuff
dndscv_input = cancereffectsizeR:::dndscv_preprocess(cesa = luad, covariates = "lung")
saveRDS(dndscv_input, "dndscv_input_single.rds")
dndscv_raw_output = lapply(dndscv_input, function(x) do.call(dndscv::dndscv, x))
# a few attributes are huge (>1 GB); drop these
dndscv_raw_output = lapply(dndscv_raw_output, function(x) { x$nbreg$terms = NULL; x$nbreg$model = NULL; x$poissmodel = NULL; return(x)})
saveRDS(dndscv_raw_output, "dndscv_raw_output_single.rds")
dndscv_out = dndscv_postprocess(cesa = luad, dndscv_raw_output = dndscv_raw_output)
saveRDS(dndscv_out@dndscv_out_list[[1]], "sel_cv.rds")
saveRDS(dndscv_out@mutrates, "mutrates.rds")

# Now run gene_mutation_rates normally
luad = gene_mutation_rates(luad, covariates = "lung")
saveRDS(luad, "cesa_for_snv.rds")

# save results to serve as expected test output
test_genes = c("EGFR", "ASXL3", "KRAS", "RYR2", "USH2A", "CSMD3", "TP53", "CSMD1", "LRP1B", 
               "ZFHX4", "FAT3", "CNTNAP5", "PCDH15", "NEB", "RYR3", "DMD", "KATNAL1", 
               "OR13H1", "KSR1")
luad = ces_snv(luad, genes = test_genes, include_nonrecurrent_variants = T)
saveRDS(luad@selection_results, "single_stage_snv_results.rds")
luad = ces_gene_epistasis(luad, genes = c("EGFR", "KRAS", "TP53"), return_all_opm_output = T)
saveRDS(luad@gene_epistasis_results, "epistasis_results.rds")
saveRDS(luad@advanced$opm_output, "epistasis_opm.rds")



anno_out = annotate_variants(dndscv_out)

saveRDS(anno_out@mutations, "mutations_anno.rds")
saveRDS(anno_out@maf, "annotated_maf_df.rds")

setwd(prev_dir)




