prev_dir = setwd(system.file("data-raw", package = "cancereffectsizeR"))

source("build_codon_point_mutation_dict.R")
codon_point_mutation_dict = build_codon_point_mutation_dict()

source("build_deconstructSigs_trinuc_string.R")
deconstructSigs_trinuc_string = build_deconstructSigs_trinuc_string()

usethis::use_data(codon_point_mutation_dict, deconstructSigs_trinuc_string, internal = TRUE, overwrite = TRUE)
setwd(prev_dir)
