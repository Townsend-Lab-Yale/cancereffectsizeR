prev_dir = setwd(system.file("data-raw", package = "cancereffectsizeR"))

source("build_deconstructSigs_stuff.R")
deconstructSigs_trinuc_string = build_deconstructSigs_trinuc_string()
deconstructSigs_notations = build_dS_notation_table()

source("build_codon_snvs_to_aa.R")
codon_sbs_to_aa = build_codon_snvs_to_aa()

source("build_codon_dbs_to_aa.R")
codon_dbs_to_aa = build_codon_dbs_to_aa()

source("build_two_codon_dbs_key.R")
two_codon_dbs_key = build_two_codon_dbs_key()


cosmic_sbs_signature_etiology = data.table::fread(system.file('extdata/cosmic_sbs_signature_summary.txt', package = 'cancereffectsizeR'))

# Pulled from COSMIC DBS mutational signature definitions
cosmic_dbs_classes = c(
  "AC>CA",
  "AC>CG",
  "AC>CT",
  "AC>GA",
  "AC>GG",
  "AC>GT",
  "AC>TA",
  "AC>TG",
  "AC>TT",
  "AT>CA",
  "AT>CC",
  "AT>CG",
  "AT>GA",
  "AT>GC",
  "AT>TA",
  "CC>AA",
  "CC>AG",
  "CC>AT",
  "CC>GA",
  "CC>GG",
  "CC>GT",
  "CC>TA",
  "CC>TG",
  "CC>TT",
  "CG>AT",
  "CG>GC",
  "CG>GT",
  "CG>TA",
  "CG>TC",
  "CG>TT",
  "CT>AA",
  "CT>AC",
  "CT>AG",
  "CT>GA",
  "CT>GC",
  "CT>GG",
  "CT>TA",
  "CT>TC",
  "CT>TG",
  "GC>AA",
  "GC>AG",
  "GC>AT",
  "GC>CA",
  "GC>CG",
  "GC>TA",
  "TA>AT",
  "TA>CG",
  "TA>CT",
  "TA>GC",
  "TA>GG",
  "TA>GT",
  "TC>AA",
  "TC>AG",
  "TC>AT",
  "TC>CA",
  "TC>CG",
  "TC>CT",
  "TC>GA",
  "TC>GG",
  "TC>GT",
  "TG>AA",
  "TG>AC",
  "TG>AT",
  "TG>CA",
  "TG>CC",
  "TG>CT",
  "TG>GA",
  "TG>GC",
  "TG>GT",
  "TT>AA",
  "TT>AC",
  "TT>AG",
  "TT>CA",
  "TT>CC",
  "TT>CG",
  "TT>GA",
  "TT>GC",
  "TT>GG"
)

cna_size_bins = data.table(pretty_label = c('L < 100 kb', '100 kb < L < 1 Mb', '1 Mb < L < 10 Mb', '10 Mb < L < 40 Mb', '40 Mb < L'),
                           cosmic_label = c('0-100kb', '100kb-1Mb', '1Mb-10Mb', '10Mb-40Mb', '>40Mb'),
                           bin_id = c('r1', 'r2', 'r3', 'r4', 'r5'),
                           bin_order = 1:5)


usethis::use_data(deconstructSigs_trinuc_string, 
                  deconstructSigs_notations, codon_sbs_to_aa, codon_dbs_to_aa,
                  two_codon_dbs_key,
                  cosmic_dbs_classes, cosmic_sbs_signature_etiology,
                  cna_size_bins,
                  internal = TRUE, overwrite = TRUE)
setwd(prev_dir)
