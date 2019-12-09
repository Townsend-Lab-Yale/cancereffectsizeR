# get_test_file and get_test_data are loaded automatically from helpers.R by testthat

test_that("ces_buildref produces expected data", {
  cds_small = get_test_data("ensembl_cds_hg19_small.rds")
  genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  refcds_and_gr_genes = suppressWarnings(ces_buildref(cds_small, genome))
  ak = get_test_data("refcds_hg19_small.rds")
  expect_equal(refcds_and_gr_genes, ak)
})

