# get_test_file and get_test_data are loaded automatically from helpers.R by testthat
tiny_maf = fread(get_test_file("tiny.hg38.maf.txt"))
test_that("load_maf and variant annotation", {
  preloaded = preload_maf(maf = tiny_maf, refset = 'ces.refset.hg38')
  expect_equivalent(preloaded, fread(get_test_file('tiny_preloaded.txt')))
  tiny = load_maf(cesa = CESAnalysis(refset = "ces.refset.hg38"), maf = preloaded)
  tiny_ak = load_cesa(get_test_file("tiny_hg38_maf_loaded.rds"))
  expect_equal(tiny$maf[order(variant_id)], tiny_ak$maf[order(variant_id)])
  expect_equal(tiny@excluded, tiny_ak@excluded, ignore.row.order = T)
  suppressWarnings({tiny_ak@mutations$amino_acid_change$constituent_sbss = NULL})
  expect_equal(tiny@mutations[c("amino_acid_change", "sbs", "aac_sbs_key")], 
               tiny_ak@mutations[c("amino_acid_change", "sbs", "aac_sbs_key")])
  expect_equal(tiny@mutations$sbs[, .N], 389)
  expect_equal(tiny@mutations$amino_acid_change[, .N], 381)
  expect_equal(tiny$samples[, .N], 9) 
   
  # same ranges should be in each coverage GenomicRange (depending on BSgenome version, little contigs may vary)
  expect_equal(lapply(tiny@coverage$exome, IRanges::ranges), lapply(tiny_ak@coverage$exome, IRanges::ranges))
  
  # Verify that calling internal function annotate_variants works the same called directly
  # Note that DBS/MNV correction occurs in load_maf, but the SBS/AACs annotated will stay the same
  variants_to_check = copy(tiny@maf)
  variants_to_check[, variant_type := NULL] # cause variants to be re-identified
  re_anno = annotate_variants(ces.refset.hg38, variants_to_check)
  expect_equal(tiny@mutations$amino_acid_change, re_anno$amino_acid_change)
  
  # column order is different when running annotate_variants() directly.
  setcolorder(re_anno$sbs, c("sbs_id", "chr", "pos", "ref", "alt", "genes", "intergenic", 
                             "trinuc_mut", "essential_splice", "nearest_pid"))
  expect_equivalent(tiny@mutations$sbs, re_anno$sbs)
  
  # expect error when adding a variant already present
  expect_error(add_variants(target_cesa = tiny, sbs_id = "10:87933130_G>C"), "all of them are already annotated")
  
  # add a variant that isn't covered by any covered regions
  tiny = add_variants(target_cesa = tiny, sbs_id = "12:132824581_A>C")
  
  
  selected = select_variants(tiny, variant_ids  = "12:132824581_A>C", genes = "TAF1C")
  expect_equal(tiny@mutations$variants_to_cov$`12:132824581_A>C`, character())
  expect_equal(selected[, .N], 2)
  selected = select_variants(tiny, min_freq = 1)
  expect_equal(sum(selected$maf_prevalence), 266) 
  expect_equal(variant_counts(tiny, "12:132824581_A>C")$N, 0)
  expect_equal(sum(variant_counts(tiny, selected$variant_id)$N), 266)
  
  # Check an essential splice site manually added to ces.refset.hg38
  tiny = add_variants(target_cesa = tiny, aac_id = 'TP53_T125T_ENSP00000269305.4')
  expect_equal(tiny$variants[variant_id == "TP53_T125T_ENSP00000269305.4", essential_splice], T)
  
  # Error when any requested gene is not in RefCDS data
  expect_error(select_variants(tiny, genes = c("KRAS", "TP53", "notagene")),
               "Some of the selected genes do not appear in the CESAnalysis reference data")
  
  # TTN is not in the data set
  expect_null(select_variants(tiny, genes = c("TTN")))
    
  
  
  # test adding covered regions and covered_regions_padding
  tiny = add_covered_regions(target_cesa = tiny, covered_regions = GRanges("chr12:132824580"), 
                      covered_regions_name = "precise_target_1", coverage_type = "targeted")
  expect_equal(tiny@mutations$variants_to_cov$`12:132824581_A>C`, character())
  
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
  expect_equal(tiny@mutations$variants_to_cov$`12:132824581_A>C`, "precise_target_2")
  
  
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


test_that("add_variants", {
  cesa = CESAnalysis('ces.refset.hg38')
  cesa = add_variants(cesa, aac_id = 'BRAF_V600E')
  expect_equal(cesa$variants[, .N], 1)
  cesa = CESAnalysis('ces.refset.hg38')
  cesa = add_variants(cesa, sbs_id = "11:18752521_G>C")
  expect_equal(cesa$variants[, .N], 1)
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

# Position 1:10001 has ambiguous trinuc context (first 10,000 bases of chr1 are N).
# FYI, the same holds for hg38.
test_that('Handle all MAF records excluded due to ambiguous trinuc context', {
  maf = as.data.table(list(Chromosome = 1, Start_Position = 10001, Reference_Allele = 'T', Tumor_Allele = 'C',
                           Unique_Patient_Identifier = 'hello'))
  expect_warning(expect_message({cesa = load_maf(cesa = CESAnalysis('ces.refset.hg19'), maf = maf)}, 
                 'excluded due to ambiguous trinucleotide'), 'Could this be whole-genome')
  expect_equal(cesa$excluded[, .N], 1)
  expect_equal(cesa$samples[, .N], 1)
  expect_equal(cesa@maf[, .N], 0)
})


