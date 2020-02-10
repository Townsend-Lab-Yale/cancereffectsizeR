test_that("multi-stage dNdScv and annotation", {
  cesa = get_test_data("cesa_for_multi_dndscv.rds")
  cesa = expect_warning(gene_level_mutation_rates(cesa, covariate_file = "lung_pca"))
  
  sel_cv = lapply(cesa@dndscv_out_list, function(x) data.table(x$sel_cv))
  sel_cv_ak = get_test_data("sel_cv_multi.rds")
  sel_cv_ak = lapply(sel_cv_ak, data.table)
  expect_equal(sel_cv, sel_cv_ak)
  
  mutrates_ak = get_test_data("mutrates_multi.rds")
  expect_equal(cesa@mutrates_list, mutrates_ak)
  
  cesa = annotate_gene_maf(cesa)
  annotated_maf = get_test_data("multi_annotated_maf_df.rds")
  expect_equal(cesa@annotated.snv.maf, annotated_maf)
})