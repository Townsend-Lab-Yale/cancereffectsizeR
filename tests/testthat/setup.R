# turn off pbapply progress bar (and save previous settings to restore on teardown)
pbopts = pbapply::pboptions(type="none")

# this overrides efforts by devtools::test to use withr::with_collate("C",...
Sys.setlocale(locale = "")

# override load_cesa warning messages
load_cesa = function(...) {
  withCallingHandlers(cancereffectsizeR::load_cesa(...), 
                      warning = function(w) {
                        if(startsWith(conditionMessage(w), "Version change")) {
                          invokeRestart("muffleWarning")
                        }
                      })
}
