get_test_data = function(filename) {
  file = system.file(paste0("tests/test_data/", filename), package = "cancereffectsizeR")
  
  # for compatibility with R 3.5, ignore the (harmless, so far) warning that comes with reading in newer RDS files
  withCallingHandlers(
    {
      object = readRDS(file)
    }, warning = function(w) {
      if (startsWith(conditionMessage(w), "cannot unserialize ALTVEC object of class 'wrap_logical'")) {
        invokeRestart("muffleWarning")
      }
    }
  )
  return(object)
}

get_test_file = function(filename) {
  path = system.file(paste0("tests/test_data/", filename), package = "cancereffectsizeR")
  return(path)
}