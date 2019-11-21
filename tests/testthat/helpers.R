get_test_data = function(filename) {
  file = system.file(paste0("tests/testthat/testdata/", filename), package = "cancereffectsizeR")
  return(readRDS(file))
}

get_test_file = function(filename) {
  path = system.file(paste0("tests/testthat/testdata/", filename), package = "cancereffectsizeR")
  return(path)
}