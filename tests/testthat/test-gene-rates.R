test_that("gene mutation rates", {
  dndscv_cesa = load_cesa(get_test_file("annotated_fruit_cesa.rds"))
  dndscv_input = get_test_data("dndscv_input_single.rds")
  dndscv_raw_output = get_test_data("dndscv_raw_output_single.rds")
  
  # Run gene_mutation_rates, but skip running dNdScv. Instead, verify that input to
  # run_dndscv matches what's expected, then load pre-generated output and continue
  dndscv_cesa = mockr::with_mock(
    run_dndscv = function(mutations, gene_list, cv, gr_genes, refdb) {
      expect_equal(list(mutations = mutations, gene_list = gene_list, cv = cv), dndscv_input)
      return(dndscv_raw_output)
    },
    gene_mutation_rates(dndscv_cesa, covariates = "lung"))
  rates_ak = get_test_data("mutrates.rds")
  expect_equal(dndscv_cesa@mutrates, rates_ak)
})