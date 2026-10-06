

get_test_data = function(filename) {
  file = testthat::test_path("..", "test_data", filename)
  
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
  path = testthat::test_path("..", "test_data", filename)
  return(path)
}