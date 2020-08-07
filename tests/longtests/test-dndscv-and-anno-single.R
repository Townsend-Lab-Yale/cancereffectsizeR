test_that("dNdScv and MAF annotation", {
  cesa = load_cesa(get_test_file("cesa_for_dndscv_and_anno.rds"))
  cesa = expect_warning(gene_mutation_rates(cesa, covariates = "lung"), "Same mutations observed in different sample")
  sel_cv = get_test_data("sel_cv.rds")
  #expect_equal(cesa@dndscv_out_list[[1]][order(sel_cv$gene_name)], sel_cv)
  mutrates = get_test_data("mutrates.rds")
  expect_equal(cesa@mutrates, mutrates)
})