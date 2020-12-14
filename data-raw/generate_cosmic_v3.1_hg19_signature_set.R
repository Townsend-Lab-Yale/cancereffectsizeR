library(data.table)
library(cancereffectsizeR)

# downloaded 08-09-20
cosmic = fread(system.file("extdata/COSMIC_SBS_v3-1.txt", package = "cancereffectsizeR"))
metadata = fread(system.file("extdata/COSMIC_v3.1_signature_metadata.txt", package = "cancereffectsizeR"))

# column names will be deconstructSigs-style trinuc mutations
dS_muts = cosmic[, deconstructSigs_notations[.(Subtype, substr(Type, 3, 3)), deconstructSigs_ID]]


# drop non-signature columns
cosmic = cosmic[, .SD, .SDcols = patterns("SBS")]
sig_names = colnames(cosmic)
cosmic[cosmic == '-'] = '0' # some empty values (reflecting something about COSMIC methods?)

# convert all to numeric (those with "-" ended up as character)
cosmic[, (sig_names) := lapply(.SD, as.numeric), .SDcols = sig_names]

cosmic_df = as.data.frame(t(cosmic))
rownames(cosmic_df) = sig_names
colnames(cosmic_df) = dS_muts

# put columns in canonical order (the order used by deconstructSigs, to avoid mistakes later)
deconstructSigs_trinuc_string = getFromNamespace("deconstructSigs_trinuc_string", "cancereffectsizeR")
cosmic_df = cosmic_df[, deconstructSigs_trinuc_string]
signature_set = list(name = "COSMIC v3.1", signatures = cosmic_df, meta = metadata)

# trigger an error if this signature set isn't valid
validate_signature_set(signature_set)

# save in hg19 reference data collection
out_path = paste0(system.file("ref_sets/ces_hg19_v1/signatures", package = "cancereffectsizeR"), '/COSMIC_v3.1_signatures.rds')
saveRDS(signature_set, out_path)
