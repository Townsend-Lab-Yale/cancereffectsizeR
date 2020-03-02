prev_dir = setwd(system.file("tests/test_data/", package = "cancereffectsizeR"))
maf_file = "luad.hg19.maf.txt"
cesa = load_maf(cesa = CESAnalysis(genome = "hg19"), maf = maf_file, sample_col = "sample_id")
small_maf = cesa@maf[Unique_Patient_Identifier %in% c("sample-15", "sample-102", "sample-30")]
# add in a recurrent mutation for testing purposes
small_maf = rbind(small_maf, list("sample-15", "1", 226180249, "G", "A"))

cesa = load_maf(cesa = CESAnalysis("hg19"), maf = small_maf)
saveRDS(cesa, "cesa_for_trinuc_weighting_calc.rds")
trimut = trinucleotide_mutation_weights(cesa)
saveRDS(trimut@trinucleotide_mutation_weights, "trinuc_mut_weighting.rds")
setwd(prev_dir)
