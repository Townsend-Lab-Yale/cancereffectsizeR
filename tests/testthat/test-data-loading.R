# get_test_file and get_test_data are loaded automatically from helpers.R by testthat
tiny_maf = get_test_file("tiny.hg19.maf.txt")
test_that("load_maf and variant annotation", {
  tiny = expect_warning(load_maf(cesa = CESAnalysis(ref_set = "ces.refset.hg19"), annotate = T, maf = tiny_maf,
                                 sample_col = "sample_id", tumor_allele_col = "Tumor_Seq_Allele2"),
                        "SNV records do not match the given reference genome")
  tiny_ak = load_cesa(get_test_file("tiny_hg19_maf_loaded.rds"))
  
  expect_equal(tiny$maf[order(variant_id)], tiny_ak$maf[order(variant_id)])
  expect_equal(tiny@excluded, tiny_ak@excluded, ignore.row.order = T)
  expect_equal(tiny@mutations, tiny_ak@mutations)
  expect_equal(tiny@mutations$snv[, .N], 269)
  expect_equal(tiny@mutations$amino_acid_change[, .N], 131)
   
  # same ranges should be in each coverage GenomicRange (depending on BSgenome version, little contigs may vary)
  expect_equal(lapply(tiny@coverage$exome, IRanges::ranges), lapply(tiny_ak@coverage$exome, IRanges::ranges))
  
  # undo annotations, verify annotate_variants works the same when called directly
  tiny@maf = tiny@maf[, .(Unique_Patient_Identifier, Chromosome, Start_Position, Reference_Allele, Tumor_Allele, variant_type)]
  tiny@maf[variant_type == "likely_mnv", variant_type := "snv"] # let annotate_variants redo the prediction
  tiny@mutations = list()
  tiny = annotate_variants(tiny)
  expect_equal(tiny@mutations, tiny_ak@mutations)
  expect_equal(tiny@maf[order(variant_id)], tiny_ak@maf[order(variant_id)])
  
  # expect error when adding a variant already present
  expect_error(add_variants(target_cesa = tiny, snv_id = "13:19752521_T>A"), "all of them are already annotated")
  
  # add a variant that isn't covered by any covered regions
  tiny = add_variants(target_cesa = tiny, snv_id = "12:132824581_C>A")
  
  
  selected = select_variants(tiny, variant_passlist  = "12:132824581_C>A", genes = "TTN")
  expect_equal(selected[variant_id == "12:132824581_C>A", unlist(covered_in)], NA_character_)
  expect_equal(selected[, .N], 3)
  selected = select_variants(tiny, variant_passlist = "12:132824581_C>A", genes = "TTN", min_freq = 0, include_subvariants = T)
  expect_equal(selected[, .N], 7)
  selected = select_variants(tiny, min_freq = 2)
  expect_equal(sum(selected$maf_frequency), 55)
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
})



test_that("load_maf edge cases", {
  # already annotated data, so should get an error when trying to load new data without annotating
  # will test the opposite situation later (see progression tests)
  tiny = load_cesa(get_test_file("tiny_hg19_maf_loaded.rds"))
  expect_error(load_maf(tiny, maf = tiny_maf, annotate = F, sample_col = "sample_id", tumor_allele_col = "Tumor_Seq_Allele2"),
               "already contains annotated variants")
  
  # you can't reload the same MAF
  expect_error(load_maf(tiny, maf = tiny_maf, annotate = T, sample_col = "sample_id", tumor_allele_col = "Tumor_Seq_Allele2"),
               "some sample IDs already appear in previously loaded data")
  
  # try loading empty/non-existent files and data
  cesa = CESAnalysis(ref_set = "ces.refset.hg19")
  expect_error(load_maf(cesa = cesa, maf = ""), "MAF not found")
  # data.table gives a warning, but load_maf throws its own error
  expect_error(expect_warning(load_maf(cesa = cesa, maf = get_test_file("empty.maf.txt")), "Input MAF data set is empty"))
  expect_error(load_maf(cesa = cesa, maf = data.frame()), "Input MAF data set is empty")
})

test_that("Sample group handling", {
  tiny = CESAnalysis(ref_set = "ces.refset.hg19")
  tiny_maf = get_test_file("tiny.hg19.maf.txt")
  
  # You can't supply a group_col to a CESAnalysis that is not stage-specific
  expect_error(load_maf(cesa = CESAnalysis(ref_set = "ces.refset.hg19"), annotate = F, maf = tiny_maf, sample_col = "sample_id", tumor_allele_col = "Tumor_Seq_Allele2", group_col = "nonexistent-column"),
               "This CESAnalysis does not specify sample groups")
  
  # If CESAnalysis is stage-specific, calls to load_maf must include group_col
  bad_maf = get_test_file("bad_stages_1.maf.txt")
  multistage = CESAnalysis(sample_groups = 1:2, ref_set = "ces.refset.hg19")
  expect_error(load_maf(multistage, maf = bad_maf), "You must specify group_col")
  
  # Ensuring that MAF data triggers error if progressions are out of bounds
  expect_error(load_maf(multistage, maf = bad_maf, group_col = "stage"), "The following groups were not declared")
  
  # Error if one sample has multiple stages listed in MAF data
  multistage = CESAnalysis(sample_groups = 1:4, ref_set = "ces.refset.hg19")
  expect_error(load_maf(multistage, maf = bad_maf, group_col = "stage"), "samples are associated with multiple groups")
  
  # Absence of a declared progression state in the data triggers a warning
  multistage = expect_warning(load_maf(multistage, annotate = F, maf = fread(bad_maf)[2:4], group_col = "stage"), "they weren't present in the MAF data")
  
  
  # Can't load with annotate = T if data has previously been loaded without annotating
  multistage = expect_error(load_maf(multistage, maf = tiny, annotate = T), "already contains unannotated records")
  
})


test_that("Coverage arguments", {
  tiny = CESAnalysis(ref_set = "ces.refset.hg19")
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

