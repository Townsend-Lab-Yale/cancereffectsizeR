# get_test_file and get_test_data are loaded automatically from helpers.R by testthat
test_that("MAF data loads correctly", {
  tiny_maf = get_test_file("tiny.hg19.maf.txt")
  tiny = load_maf(cesa = CESAnalysis(), maf = tiny_maf, sample_col = "sample_id", tumor_allele_col = "Tumor_Seq_Allele2")
  tiny_ak = get_test_data("tiny_hg19_maf_loaded.rds")
  expect_identical(data.frame(tiny@maf), data.frame(tiny_ak@maf))
  expect_identical(data.frame(tiny@excluded), data.frame(tiny_ak@excluded))
  expect_identical(tiny@coverage, tiny_ak@coverage)
  expect_error(load_maf(tiny, maf = tiny_maf, sample_col = "sample_id", tumor_allele_col = "Tumor_Seq_Allele2"),
               "Sample identifiers in new data have overlap")
  
  # try loading empty/non-existent file/data
  cesa = CESAnalysis()
  expect_error(load_maf(cesa = cesa, maf = ""), "MAF not found")
  # data.table gives a warning, but load_maf throws its own error
  expect_error(expect_warning(load_maf(cesa = cesa, maf = get_test_file("empty.maf.txt")), "Input MAF data set is empty"))
  expect_error(load_maf(cesa = cesa, maf = data.frame()), "Input MAF data set is empty")
})

test_that("Progression stage handling", {
  tiny = CESAnalysis()
  tiny_maf = get_test_file("tiny.hg19.maf.txt")
  
  # You can't supply a progression_col to a CESAnalysis that is not stage-specific
  expect_error(load_maf(cesa = CESAnalysis(), maf = tiny_maf, sample_col = "sample_id", tumor_allele_col = "Tumor_Seq_Allele2", progression_col = "nonexistent-column"),
               "This CESAnalysis is not stage-specific")
  
  # If CESAnalysis is stage-specific, calls to load_maf must include progression_col
  bad_maf = get_test_file("bad_stages_1.maf.txt")
  multistage = CESAnalysis(progression_order = 1:2)
  expect_error(load_maf(multistage, maf = bad_maf), "You must supply a progression_col")
  
  # Ensuring that MAF data triggers error if progressions are out of bounds
  expect_error(load_maf(multistage, maf = bad_maf, progression_col = "stage"), "Unexpected progressions stage")
  
  # Error if one sample has multiple stages listed in MAF data
  multistage = CESAnalysis(progression_order = 1:4)
  expect_error(load_maf(multistage, maf = bad_maf, progression_col = "stage"), "multiple progression stages")
})
