# get_test_file and get_test_data are loaded automatically from helpers.R by testthat
tiny_maf = fread(get_test_file("tiny.hg19.maf.txt"))
test_that("load_maf and variant annotation", {
  tiny = load_maf(cesa = CESAnalysis(refset = "ces.refset.hg19"), maf = tiny_maf)
  tiny_ak = load_cesa(get_test_file("tiny_hg19_maf_loaded.rds"))
  expect_equal(tiny$maf[order(variant_id)], tiny_ak$maf[order(variant_id)])
  expect_equal(tiny@excluded, tiny_ak@excluded, ignore.row.order = T)
  expect_equal(tiny@mutations, tiny_ak@mutations)
  expect_equal(tiny@mutations$snv[, .N], 269)
  expect_equal(tiny@mutations$amino_acid_change[, .N], 131)
  expect_equal(tiny$samples[, .N],84) 
   
  # same ranges should be in each coverage GenomicRange (depending on BSgenome version, little contigs may vary)
  expect_equal(lapply(tiny@coverage$exome, IRanges::ranges), lapply(tiny_ak@coverage$exome, IRanges::ranges))
  
  # Verify that calling internal function annotate_variants works the same called directly
  # Note that DBS/MNV correction occurs in load_maf, but the SNVs/AACs annotated will stay the same
  variants_to_check = copy(tiny@maf)
  variants_to_check[, variant_type := NULL] # cause variants to be re-identified
  re_anno = annotate_variants(ces.refset.hg19, variants_to_check)
  expect_equal(tiny@mutations$amino_acid_change[, -"covered_in"], re_anno$amino_acid_change)
  expect_equal(tiny@mutations$snv[, -"covered_in"], re_anno$snv)
  
  # expect error when adding a variant already present
  expect_error(add_variants(target_cesa = tiny, snv_id = "13:19752521_T>A"), "all of them are already annotated")
  
  # add a variant that isn't covered by any covered regions
  tiny = add_variants(target_cesa = tiny, snv_id = "12:132824581_C>A")
  
  
  selected = select_variants(tiny, variant_ids  = "12:132824581_C>A", genes = "TTN")
  expect_equal(selected[variant_id == "12:132824581_C>A", unlist(covered_in)], NA_character_)
  expect_equal(selected[, .N], 3)
  selected = select_variants(tiny, variant_ids = "12:132824581_C>A", genes = "TTN", min_freq = 0, include_subvariants = T)
  expect_equal(selected[, .N], 7)
  selected = select_variants(tiny, min_freq = 2)
  expect_equal(sum(selected$maf_prevalence), 55) # if you get 57, MNV on sample-D probably got included
  expect_equal(variant_counts(tiny, "12:132824581_C>A")$total_prevalence, 0)
  expect_equal(sum(variant_counts(tiny, selected$variant_id)$total_prevalence), 55)
  expect_equal(selected[variant_id == "TP53_T125T_ENSP00000269305", essential_splice], T)
  
  # Error when any requested gene is not in RefCDS data
  expect_error(select_variants(tiny, genes = c("KRAS", "TP53", "notagene")),
               "Some of the selected genes do not appear in the CESAnalysis reference data")
  
  # AC006486.1 is not in the data set
  expect_null(select_variants(tiny, genes = c("AC006486.1")))
    
  
  
  # test adding covered regions and covered_regions_padding
  tiny = add_covered_regions(target_cesa = tiny, covered_regions = GRanges("chr12:132824580"), 
                      covered_regions_name = "precise_target_1", coverage_type = "targeted")
  expect_equal(tiny@mutations$snv["12:132824581_C>A", unlist(covered_in)], character())
  
  # load again, with different ranges
  expect_error(add_covered_regions(target_cesa = tiny, covered_regions = GRanges("chr12:100"), 
                             covered_regions_name = "precise_target_1", coverage_type = "targeted"),
                             "covered_regions do not exactly match")
  
  # load again, same ranges, but call it "exome"
  expect_error(add_covered_regions(target_cesa = tiny, covered_regions = GRanges("chr12:132824580"), 
                      covered_regions_name = "precise_target_1", coverage_type = "exome"), "already been used")
  
  # bad covered_regions_name
  expect_error(add_covered_regions(target_cesa = tiny, covered_regions = GRanges("chr12:1-100"),
                      covered_regions_name = "bad***", coverage_type = "exome"), "Invalid covered_regions_name")
  
  # load again with padding
  tiny = add_covered_regions(target_cesa = tiny, covered_regions = GRanges("chr12:132824580"), 
                      covered_regions_name = "precise_target_2", coverage_type = "targeted", covered_regions_padding = 10)
  expect_equal(tiny@mutations$snv["12:132824581_C>A", unlist(covered_in)], "precise_target_2")
  
  
  # check load_sample_data and clear_sample_data
  original_samples = copy(tiny@samples)
  t1 = original_samples[, .(Unique_Patient_Identifier = rev(Unique_Patient_Identifier), 
                            coverage, covered_regions, num = 1:.N)]
  new_info = load_sample_data(tiny, t1)@samples
  
  expect_equal(new_info[, -"num"], original_samples, check.attributes = F)
  
  altered = copy(t1)
  altered$coverage[5] = 'hello'
  expect_error(load_sample_data(tiny, altered), "data don't match")
  
  t2 = t1[1:5]
  new_info = load_sample_data(tiny, t2)@samples
  expect_identical(new_info$num, c(rep(NA, new_info[, .N - 5]), 5:1))
  
  expect_error(load_sample_data(tiny, t1[, -"Unique_Patient_Identifier"]), "must have a Unique_Patient_Identifier")
  expect_error(load_sample_data(tiny, t1[, .(Unique_Patient_Identifier, coverage)]), "There are no new data columns")
  
  new_info[, num2 := 1:.N]
  tiny = load_sample_data(tiny, new_info)
  
  expect_error(load_sample_data(tiny, new_info), 'non-missing values overwritten by sample_data')
  expect_error(clear_sample_data(tiny, c('num', 'num2', 'Unique_Patient_Identifier')), 'Internal columns')
  tiny = clear_sample_data(tiny, c('num', 'num2'))
  all.equal(tiny@samples, original_samples)
})



test_that("load_maf edge cases", {
  tiny = load_cesa(get_test_file("tiny_hg19_maf_loaded.rds"))
  
  # you can't reload the same MAF (suppressing reference allele mismatch warning)
  expect_error(expect_warning(load_maf(tiny, maf = tiny_maf)),
               "some sample IDs already appear in previously loaded data")
  
  # try loading empty/non-existent files and data
  cesa = CESAnalysis(refset = "ces.refset.hg19")
  expect_error(load_maf(cesa = cesa, maf = ""), "MAF not found")
  # data.table gives a warning, but load_maf throws its own error
  expect_error(expect_warning(load_maf(cesa = cesa, maf = get_test_file("empty.maf.txt")), "Input MAF data set is empty"))
  expect_error(load_maf(cesa = cesa, maf = data.frame()), "Input MAF data set is empty")
})

test_that("Sample group handling", {
  tiny = CESAnalysis(refset = ces.refset.hg19) # note use of ces.refset.hg19 directly (also allowed)
  tiny_maf = get_test_file("tiny.hg19.maf.txt")
  
  # You can't supply a group_col to a CESAnalysis that does not have predefined sample groups
  # Also, group_col is depcrecasted and will be removed soon.
  expect_error(load_maf(cesa = CESAnalysis(refset = "ces.refset.hg19"), maf = tiny_maf, group_col = "nonexistent-column"),
               "group_col is deprecated")
  
  # If CESAnalysis is stage-specific, calls to load_maf must include group_col
  bad_maf = get_test_file("bad_stages_1.maf.txt")
  multistage = expect_warning(CESAnalysis(sample_groups = 1:2, refset = "ces.refset.hg19"), 'sample_groups is deprecated')
  
  # Suppress deprecation warning
  expect_error(suppressWarnings(load_maf(multistage, maf = bad_maf)), "You must specify group_col")
  
  # Suppress deprecation warning
  # Ensuring that MAF data triggers error if progressions are out of bounds
  expect_error(suppressWarnings(load_maf(multistage, maf = bad_maf, group_col = "stage")), "The following groups were not declared")
  
  # Error if one sample has multiple stages listed in MAF data
  multistage = expect_warning(CESAnalysis(sample_groups = 1:4, refset = "ces.refset.hg19"), 'sample_groups is deprecated')
  expect_error(load_maf(multistage, maf = bad_maf, group_col = "stage"), "samples are associated with multiple groups")
  
  # Absence of a declared progression state in the data triggers a notification
  # Also, this test happens to make sure annotate_variants doesn't fail when there are zero possible splice site records
  multistage = expect_message(suppressWarnings(load_maf(multistage, maf = fread(bad_maf)[2:4], group_col = "stage")), "they weren't present in the MAF data")
})


test_that("Coverage arguments", {
  tiny = CESAnalysis(refset = "ces.refset.hg19")
  maf = data.table() # coverage arguments get validated before data actually loaded in load_maf 
  expect_error(load_maf(tiny, maf = maf, coverage = c("genome", "hi")), 
               "coverage must be \"exome")
  expect_error(load_maf(tiny, maf = maf, coverage = NULL),
               "coverage must be \"exome")
  expect_error(load_maf(tiny, maf = maf, covered_regions_name = "exome"),
               "covered_regions_name was supplied, but covered_regions wasn't")
  expect_error(load_maf(tiny, maf = maf, covered_regions = data.table()),
               "You must supply a name for your covered_regions")
  expect_error(load_maf(tiny, maf = maf, coverage = "targeted"), 
               "can't load targeted data without covered_regions")
})

