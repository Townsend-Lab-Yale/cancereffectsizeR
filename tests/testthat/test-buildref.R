# get_test_file and get_test_data are loaded automatically from helpers.R by testthat

test_that("ces_buildref produces expected data", {
  cds_small = fread(get_test_file("cds_small_stop_codons_included.txt"))
  refcds_and_gr_genes = build_RefCDS(cds_small, genome = "hg19", cds_ranges_lack_stop_codons = F)
  ak = get_test_data("refcds_hg19_small.rds")
  expect_equal(refcds_and_gr_genes[[1]], ak[[1]]) # RefCDS (currently broken on switch to R 4.0)
  expect_equal(refcds_and_gr_genes[[2]], ak[[2]], check.attributes = F) # gr_genes 
})

