prev_dir = setwd(system.file("tests/test_data/", package = "cancereffectsizeR"))
maf_file = "luad.hg19.maf.txt"
cesa = load_maf(cesa = CESAnalysis("ces_hg19_v1"), maf = maf_file, sample_col = "sample_id", annotate = F)

cesa = trinuc_mutation_rates(cesa, cores = 4, signature_set = "COSMIC_v3.1",
                             signatures_to_remove = suggest_cosmic_signatures_to_remove("LUAD", treatment_naive = TRUE, quiet = TRUE))

fwrite(cesa$trinuc_rates, "luad_hg19_trinuc_rates.txt", sep = "\t")
fwrite(cesa$mutational_signatures, "luad_hg19_sig_table.txt", sep = "\t")

## Create small trinuc weight data for tests
# These samples are a mix of LUAD and fabricated data; one tumor is special because dS returns 100% artifact weightings on it
# There is one recurrent SNV that appears in both good_A and good_B
cesa = load_maf(cesa = CESAnalysis("ces_hg19_v1"), maf = get_test_file("trinuc_rate_test_maf.txt"), annotate = F)
saveRDS(cesa, "cesa_for_trinuc_weighting_calc.rds")
trimut = trinuc_mutation_rates(cesa, signature_set = "COSMIC_v3.1", 
                               signatures_to_remove = suggest_cosmic_signatures_to_remove("LUAD", TRUE, TRUE))
saveRDS(trimut@trinucleotide_mutation_weights, "trinuc_mut_weighting.rds")

