# turn off pbapply progress bar (and save previous settings to restore on teardown)
pbopts = pbapply::pboptions(type="none")

# restore user preference on teardown
teardown(pbapply::pboptions(pbopts))

# this overrides efforts by devtools::test to use withr::with_collate("C",...
Sys.setlocale(locale = "")

# maximize text width so that messages don't get newlines inserted,
# which otherwise makes expect_message fail sometimes
width_opt = options(width = 10000)$width

# restore text width on teardown
teardown(options(width = width_opt))

# override load_cesa warning messages
load_cesa = function(...) {
  withCallingHandlers(cancereffectsizeR::load_cesa(...), 
                      warning = function(w) {
                        if(startsWith(conditionMessage(w), "Version change")) {
                          invokeRestart("muffleWarning")
                        }
                      })
}
