# turn off pbapply progress bar (and save previous settings to restore on teardown)
pbopts = pbapply::pboptions(type="none")

# this overrides efforts by devtools::test to use withr::with_collate("C",...
Sys.setlocale(locale = "")
