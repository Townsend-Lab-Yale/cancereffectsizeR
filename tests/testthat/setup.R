# turn off pbapply progress bar (and save previous settings to restore on teardown)
pbopts = pbapply::pboptions(type="none")

# restore user preference on teardown
withr::defer(pbapply::pboptions(pbopts), teardown_env())

# this overrides efforts by devtools::test to use withr::with_collate("C",...
Sys.setlocale(locale = "")

# Set text width so that messages don't get unexpected newlines inserted,
# which can expect_message fail
width_opt = options(width = 150)$width

# restore text width on teardown
withr::defer(options(width = width_opt), teardown_env())


# override load_cesa warning messages
load_cesa = function(...) {
  withCallingHandlers(cancereffectsizeR::load_cesa(...), 
                      warning = function(w) {
                        if(startsWith(conditionMessage(w), "Version change")) {
                          invokeRestart("muffleWarning")
                        }
                      })
}
