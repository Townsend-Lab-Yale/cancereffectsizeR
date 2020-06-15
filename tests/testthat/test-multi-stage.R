# note that mutational signature calculation doesn't use stage information, so no tests for that step here
# testing the pre- and post-dndscv processing functions, but not actually running dndscv (only run in long tests)

test_that("multi-stage dNdScv and annotation", {
  cesa = get_test_data("cesa_for_multi_dndscv.rds")
  dndscv_input = cancereffectsizeR:::dndscv_preprocess(cesa = cesa, covariate_file = "lung_pca")
  dndscv_input_ak = get_test_data("dndscv_input_multi.rds")
  # leaving out the refdb temporary filename, which will vary
  expect_equal(lapply(dndscv_input, function(x) x[-4]), lapply(dndscv_input_ak, function(x) x[-4])) # refdb path will very on dev/prod due to inst dir
  
  dndscv_output = get_test_data("dndscv_raw_output_multi.rds")
  cesa = dndscv_postprocess(cesa, dndscv_output)
  
  sel_cv = lapply(cesa@dndscv_out_list, function(x) data.table(x$sel_cv))
  sel_cv_ak = get_test_data("sel_cv_multi.rds")
  sel_cv_ak = lapply(sel_cv_ak, data.table)
  expect_equal(sel_cv, sel_cv_ak)
  
  mutrates_ak = get_test_data("mutrates_multi.rds")
  expect_equal(cesa@mutrates, mutrates_ak)
  
  ## Annotation currently ignores stage, so single-stage testing sufficient
  # cesa = annotate_variants(cesa)
  # annotated_maf = get_test_data("multi_annotated_maf_df.rds")
  # expect_equal(cesa@maf, annotated_maf)
})



test_that("multi-stage SNV effect size calculation", {
  cesa = get_test_data("cesa_for_snv_multi.rds")
  test_genes = c("TTN", "KRAS", "RYR2", "EGFR", "TP53", "ASXL3","IFITM2")
  cesa = ces_snv(cesa, genes = test_genes, include_nonrecurrent_variants = T)
  results = cesa@selection_results
  results_ak = get_test_data("multi_stage_snv_results.rds")
  expect_equal(results_ak, results, tolerance = 1e-5)
})
