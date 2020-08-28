# get_test_file and get_test_data are loaded automatically from helpers.R by testthat
tiny_maf = get_test_file("tiny.hg19.maf.txt")
test_that("load_maf and variant annotation", {
  tiny = expect_warning(load_maf(cesa = CESAnalysis(genome="hg19"), annotate = T, maf = tiny_maf,
                                 sample_col = "sample_id", tumor_allele_col = "Tumor_Seq_Allele2"),
                        "SNV records do not match the given reference genome")
  tiny_ak = load_cesa(get_test_file("tiny_hg19_maf_loaded.rds"))
  
  expect_equal(tiny$maf, tiny_ak$maf)
  expect_equal(tiny@excluded, tiny_ak@excluded)
  expect_equal(tiny@mutations, tiny_ak@mutations)
  expect_equal(tiny@mutations$snv[, .N], 264)
  expect_equal(tiny@mutations$amino_acid_change[, .N], 128)
   
  # same ranges should be in each coverage GenomicRange (depending on BSgenome version, little contigs may vary)
  expect_equal(lapply(tiny@coverage$exome, IRanges::ranges), lapply(tiny_ak@coverage$exome, IRanges::ranges))
  
  # undo annotations, verify annotate_variants works the same when called directly
  tiny@maf = tiny@maf[, .(Unique_Patient_Identifier, Chromosome, Start_Position, Reference_Allele, Tumor_Allele, Variant_Type)]
  tiny@mutations = list()
  tiny = annotate_variants(tiny)
  expect_equal(tiny@mutations, tiny_ak@mutations)
  expect_equal(tiny@maf, tiny_ak@maf)
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
  cesa = CESAnalysis(genome = "hg19")
  expect_error(load_maf(cesa = cesa, maf = ""), "MAF not found")
  # data.table gives a warning, but load_maf throws its own error
  expect_error(expect_warning(load_maf(cesa = cesa, maf = get_test_file("empty.maf.txt")), "Input MAF data set is empty"))
  expect_error(load_maf(cesa = cesa, maf = data.frame()), "Input MAF data set is empty")
})

test_that("Progression stage handling", {
  tiny = CESAnalysis(genome="hg19")
  tiny_maf = get_test_file("tiny.hg19.maf.txt")
  
  # You can't supply a progression_col to a CESAnalysis that is not stage-specific
  expect_error(load_maf(cesa = CESAnalysis(genome="hg19"), annotate = F, maf = tiny_maf, sample_col = "sample_id", tumor_allele_col = "Tumor_Seq_Allele2", progression_col = "nonexistent-column"),
               "This CESAnalysis does not incorporate tumor progression")
  
  # If CESAnalysis is stage-specific, calls to load_maf must include progression_col
  bad_maf = get_test_file("bad_stages_1.maf.txt")
  multistage = CESAnalysis(progression_order = 1:2, genome = "hg19")
  expect_error(load_maf(multistage, maf = bad_maf), "You must specify progression_col")
  
  # Ensuring that MAF data triggers error if progressions are out of bounds
  expect_error(load_maf(multistage, maf = bad_maf, progression_col = "stage"), "The following progressions were not declared")
  
  # Error if one sample has multiple stages listed in MAF data
  multistage = CESAnalysis(progression_order = 1:4, genome = "hg19")
  expect_error(load_maf(multistage, maf = bad_maf, progression_col = "stage"), "samples are associated with multiple progressions")
  
  # Absence of a declared progression state in the data triggers a warning
  multistage = expect_warning(load_maf(multistage, annotate = F, maf = fread(bad_maf)[2:4], progression_col = "stage"), "they weren't present in the MAF data")
  
  # Can't load with annotate = T if data has previously been loaded without annotating
  multistage = expect_error(load_maf(multistage, maf = tiny, annotate = T), "already contains unannotated records")
})


test_that("Coverage arguments", {
  tiny = CESAnalysis("hg19")
  maf = data.table() # coverage arguments get validated before data actually loaded in load_maf 
  expect_error(load_maf(tiny, maf = maf, coverage = c("genome", "hi")), 
               "coverage must be \"exome")
  expect_error(load_maf(tiny, maf = maf, coverage = NULL),
               "coverage must be \"exome")
  expect_error(load_maf(tiny, maf = maf, covered_regions_name = "exome"),
               "covered_regions_name was supplied, but covered_regions wasn't")
  expect_error(load_maf(tiny, maf = maf, covered_regions = data.table()),
               "You must supply a name for your covered_regions")
  expect_error(load_maf(tiny, maf = maf, coverage = "genome", covered_regions = data.table(), covered_regions_name = "hi"),
               "covered_regions should be left NULL when coverage is \"genome")
  expect_error(load_maf(tiny, maf = maf, coverage = "targeted"), 
               "can't load targeted data without covered_regions")
})

