prev_dir = setwd(system.file("data-raw", package = "cancereffectsizeR"))

source("build_deconstructSigs_stuff.R")
deconstructSigs_trinuc_string = build_deconstructSigs_trinuc_string()
deconstructSigs_notations = build_dS_notation_table()

cosmic_v3_signature_metadata = data.table::fread("COSMIC_v3_signature_metadata.txt")

source("build_codon_snvs_to_aa.R")
codon_snvs_to_aa = build_codon_snvs_to_aa()

usethis::use_data(deconstructSigs_trinuc_string, 
                  deconstructSigs_notations, cosmic_v3_signature_metadata, codon_snvs_to_aa,
                  internal = TRUE, overwrite = TRUE)
setwd(prev_dir)
