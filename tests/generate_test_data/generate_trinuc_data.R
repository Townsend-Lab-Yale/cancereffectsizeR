prev_dir = setwd(system.file("tests/test_data/", package = "cancereffectsizeR"))
maf_file = "luad.hg19.maf.txt"
cesa = load_maf(cesa = CESAnalysis(genome = "hg19"), maf = maf_file, sample_col = "sample_id")

cesa = trinuc_mutation_rates(cesa, cores = 4, 
                             signatures_to_remove = suggest_cosmic_v3_signatures_to_remove("LUAD", treatment_naive = TRUE, quiet = TRUE))
saveRDS(cesa@trinucleotide_mutation_weights, "luad_trinucleotide_mutation_weights.rds")

## Create small trinuc weight data for tests
# these samples have been chosen because they (randomly) have a fast runtime in the group average dS step
small_maf = cesa$maf[Unique_Patient_Identifier %in% c("sample-102", "sample-1", "sample-31")]
# add in a recurrent mutation for testing purposes
small_maf = rbind(small_maf, list("sample-1", "7", 124404124, "G", "C", "SNV"))

cesa = load_maf(cesa = CESAnalysis("hg19"), maf = small_maf)
saveRDS(cesa, "cesa_for_trinuc_weighting_calc.rds")
trimut = trinuc_mutation_rates(cesa, signatures_to_remove = suggest_cosmic_v3_signatures_to_remove("LUAD", TRUE, TRUE))
saveRDS(trimut@trinucleotide_mutation_weights, "trinuc_mut_weighting.rds")


setwd(prev_dir)
