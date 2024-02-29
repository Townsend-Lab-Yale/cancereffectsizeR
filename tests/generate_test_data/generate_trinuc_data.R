prev_dir = setwd(system.file("tests/test_data/", package = "cancereffectsizeR"))

## Create small trinuc weight data for tests
# These samples are a mix of LUAD and fabricated data; one tumor is special because dS returns 100% artifact weightings on it
# There is one recurrent SNV that appears in both good_A and good_B
cesa = load_maf(cesa = CESAnalysis("ces.refset.hg19"), maf = get_test_file("trinuc_rate_test_maf.txt"))
save_cesa(cesa, "cesa_for_trinuc_weighting_calc.rds")
trimut = trinuc_mutation_rates(cesa, signature_set = "COSMIC_v3.1", 
                               signature_exclusions = suggest_cosmic_signature_exclusions("LUAD", TRUE, TRUE),
                               signature_extractor = 'deconstructSigs')
saveRDS(trimut@trinucleotide_mutation_weights, "trinuc_mut_weighting.rds")
setwd(prev_dir)

