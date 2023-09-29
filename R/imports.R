#' @import data.table
#' @import GenomeInfoDb
#' @import BSgenome
#' @importFrom utils .DollarNames data packageVersion download.file tail
#' @importFrom methods is new
#' @importFrom stats na.omit predict setNames dist window

.datatable.aware = TRUE
.ces_ref_data = new.env()
options(datatable.prettyprint.char = 40) 

# Data package names and minimum required version
.official_refsets = list(ces.refset.hg19 = as.package_version("1.1.0"), ces.refset.hg38 = as.package_version("1.2.0"))

# If refset packages are loaded, put their data in .ces_ref_data for easy access
for(refset in names(.official_refsets)) {
  if(refset %in% loadedNamespaces()) {
    .ces_ref_data[[refset]] = get(refset, envir = as.environment(paste0('package:', refset)))
  }
}

snv_annotation_template = data.table(snv_id = character(), chr = character(), pos = numeric(), 
                                     ref = character(), alt = character(), genes = list(), intergenic = logical(), 
                                     trinuc_mut = character(), essential_splice = logical(), 
                                     nearest_pid = list())

aac_annotation_template = data.table(aac_id = character(), gene = character(), aachange = character(), 
                       strand = integer(), chr = character(), pid = character(), aa_ref = character(), 
                       aa_pos = numeric(), aa_alt = character(), 
                       nt1_pos = numeric(), nt2_pos = numeric(), nt3_pos = numeric(), 
                       coding_seq = character(), essential_splice = logical())


aac_snv_key_template = data.table(aac_id = character(), snv_id = character(), multi_anno_site = logical(), key = 'aac_id')


dbs_annotation_template = data.table(dbs_id = character(), chr = character(), pos = numeric(), ref = character(), 
                                     alt = character(), intergenic = character(), essential_splice = character(),
                                     cosmic_dbs_class = character())

dbs_codon_change_template = data.table(dbs_aac_id = character(), chr = character(), pid = character(), 
                                       essential_splice = character(), strand = character(), gene = character(), 
                                       aa_ref = character(), aa_pos = character(), nt1_pos = character(), 
                                       nt2_pos = character(), nt3_pos = character(), coding_seq = character(), 
                                       aa_alt = character(), aachange = character())

aac_dbs_key_template = data.table(dbs_id = character(), dbs_aac_id = character(), multi_anno_site = logical())

sample_table_template = data.table(Unique_Patient_Identifier = character(), coverage = character(), 
                                  covered_regions = character(), sig_analysis_grp = integer(), gene_rate_grp = integer(),
                                  maf_source = character())

# for use when identifying a column previously handled by cancereffectsizeR
# "merged_into_other_variant" replaced with "merged_with_nearby_variant" in 2.6.4
preload_problems = c('missing_values', 'not_variant', 'duplicate_record', 'failed_liftOver', 
                     'duplicate_record_after_liftOver', 'unsupported_chr', 'out_of_bounds', 
                     'reference_mismatch', "merged_into_dbs_variant", "merged_into_other_variant",
                     "duplicate_from_TCGA_sample_merge", "merged_with_nearby_variant", "invalid_record")


# format a string the way R should automatically, then feed it to message()
pretty_message = function(msg, emit = T, black = emit) {
  msg = paste0(strwrap(msg), collapse = "\n")
  if (black) {
    # If current theme is dark, make the message white
    if(Sys.getenv("RSTUDIO") == "1" && rstudioapi::getThemeInfo()$dark) {
      msg = crayon::white(msg)
    } else {
      msg = crayon::black(msg)
    }
  }
  if (emit) {
    message(msg)
  } else {
    return(msg)
  }
}

NULL
