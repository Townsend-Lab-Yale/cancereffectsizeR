prev_dir = setwd(system.file("tests/test_data/", package = "cancereffectsizeR"))

## Create small trinuc weight data for tests
# There is one recurrent SNV that appears in both good_A and good_B
# The sample called zeroed-out was special in that deconstructSigs gave it all artifact weights. However,
# MutationalPatterns is now used, and it no longer gets all artifact weights. If we find another such sample,
# we will add it to testing.
cesa = load_maf(cesa = CESAnalysis("ces.refset.hg19"), maf = get_test_file("trinuc_rate_test_maf.txt"))
save_cesa(cesa, "cesa_for_trinuc_weighting_calc.rds")
trimut = trinuc_mutation_rates(cesa, signature_set = "COSMIC_v3.1", 
                               signature_exclusions = suggest_cosmic_signature_exclusions("LUAD", TRUE, TRUE),
                               signature_extractor = 'MutationalPatterns')
saveRDS(trimut@trinucleotide_mutation_weights, "trinuc_mut_weighting.rds")
setwd(prev_dir)

