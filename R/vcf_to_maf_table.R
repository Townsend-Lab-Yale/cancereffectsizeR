#' Read a VCF into an MAF-like table
#' 
#' This function loads VCF files into MAF-like tables. A Tumor_Sample_Barcode column is
#' added, and the contents of the POS/REF/ALT fields are converted to match the style used
#' by MAF files for Start_Position/Reference_Allele/Tumor_Seq_Allele2. The VCF file
#' should represent high-confidence somatic variant calls for a single tumor sample.
#' 
#' @param vcfs Vector of VCF file paths, or a list of VCF-like data.tables, or a single data.table.
#' @param sample_ids Identifiers to populate Tumor_Sample_Barcode, one per VCF.
#' @return A single data.table with MAF-style fields, suitable for use with cancereffectsizeR.
#' @export
vcfs_to_maf_table = function(vcfs, sample_ids) {
  if (is(vcfs, 'list') && length(vcfs) > 0) {
    if(! all(sapply(vcfs, is.data.frame))) {
      stop('vcfs should be a vector of file paths or a list of VCF-like data.tables (or a single data.table).')
    }
  } else if(is.character(vcfs) && length(vcfs) > 0) {
    if(uniqueN(vcfs) != length(vcfs)) {
      stop("Same file(s) specified multiple times in vcfs.")
    }
    nonexistent = vcfs[! file.exists(vcfs)]
    if(length(nonexistent) > 0) {
      msg = paste0('Some input VCF(s) do not exist: ', paste(nonexistent, collapse = ", "), '.')
      stop(pretty_message(msg, emit = F))
    }
  } else if (! is.data.frame(vcfs)) {
    stop("vcfs should be a vector of file paths or a list of VCF-like data.tables (or a single data.table).")
  } else {
    vcfs = list(vcfs) # read_vcf expects file paths or a list
  }
  if(!is.character(sample_ids)) {
    stop("sample_ids should be character vector, with one ID per VCF.")
  }
  
  legal_sample_id_pattern = '^[a-z0-9][0-9a-z\\_\\-\\.:/]*$'
  illegal_ids = sample_ids[! grepl(legal_sample_id_pattern, tolower(sample_ids), perl = T)]
  num_illegal = length(illegal_ids)
  if (num_illegal > 0) {
    if(num_illegal > 50) {
      illegal_ids = c(illegal_ids[1:40], paste0("and ", num_illegal - 40, " others"))
    }
    msg = paste0("sample_id should start with letter/number and use only letters, numbers, and the characters '-', '_', '.', ':', '/'. ",
            num_illegal, " invalid IDs: ", paste(illegal_ids, collapse = ", "), '.')
    stop(pretty_message(msg, emit = F))
  }
  
  if(length(sample_ids) != uniqueN(sample_ids)) {
    stop("Same sample ID specified more multiple VCFs.")
  }
  
  if(length(sample_ids) != length(vcfs)) {
    stop("sample_ids is not same length as vcfs.")
  }
  
  # Identify each VCF with file path (basename if possible) or list index.
  if(is.character(vcfs)) {
    vcf_names = basename(vcfs)
    if(uniqueN(vcf_names) != length(vcf_names)) {
      vcf_names = vcfs
    }
    output_names = vcf_names
  } else {
    vcf_names = paste0('(VCF ', 1:length(sample_ids), ' in list)')
    output_names = 1:length(sample_ids)
  }
  
  message("Reading ", length(vcfs), ' VCFs...')
  maf = mapply(read_vcf, vcf = vcfs, sample_id = sample_ids, 
                         vcf_name = vcf_names, SIMPLIFY = FALSE)
  names(maf) = output_names
  maf = rbindlist(maf, idcol = 'source')
  
  missing_filter_col = maf[is.na(FILTER), unique(Tumor_Sample_Barcode)]
  if(length(missing_filter_col) > 0) {
    msg = paste0("The FILTER field was missing from ", length(missing_filter_col), " input VCFs. ",
                 "FILTER has been set to NA in the output MAF table for these samples. To ",
                 "ensure that VCF-to-MAF conversion was appropriate, verify that essential fields (CHR, POS, REF, ALT) match VCF ",
                 "specifications and that the data was appropriately filtered to high-confidence calls.")
    warning(pretty_message(msg, emit =F))
  }
  
  filter_values = maf[! is.na(FILTER), unique(FILTER)]
  if(identical(filter_values, '.')) {
    msg = paste0("FILTER field is all '.' (as opposed to PASS). Verify that input VCFs consist of ",
                 "only high-confidence somatic calls. (Depending on how VCFs were made, there may or may not be a problem.)")
    warning(pretty_message(msg, emit = F))
  } else if(! identical(filter_values, 'PASS') && length(filter_values) > 0) {
      other_filter_values = setdiff(filter_values, 'PASS')
      other_filter_values = unique(unlist(strsplit(other_filter_values, ';')))
      other_filter_values[other_filter_values == '.'] = '"."'
      
      msg = paste0("Some VCF(s), and the resulting MAF output, have values besides 'PASS' in the FILTER field (", 
                   paste(other_filter_values, collapse = ', '), '). ',
                   "Make sure to only include high-confidence variant calls in downstream analyses.")
      warning(pretty_message(msg, emit = F))
  }
  
  maybe_multi_sample = maf[possible_multi_sample == T, unique(Tumor_Sample_Barcode)]
  num_multi = length(maybe_multi_sample)
  if(num_multi > 0) {
    if(num_multi > 50) {
      maybe_multi_sample = c(maybe_multi_sample[1:40], paste0('and ', num_multi - 40, ' more'))
    }
    msg = paste0(num_multi, " VCFs had extra fields: " , paste(maybe_multi_sample, collapse = ", "), '. ',
                 "If these fields represent additional tumor/normal samples ",
                 "(i.e., the files cover more than single tumor samples or tumor/normal pairs), find another way to convert your data, ",
                 "as this function does not support multi-sample VCFs.")
    warning(pretty_message(msg, emit = F))
  }
  
  with_complex = maf[has_complex == T, unique(Tumor_Sample_Barcode)]
  num_complex = length(with_complex)
  if(num_complex > 0) {
    if(num_complex > 50) {
      with_complex = c(with_complex[1:40], paste0('and ', num_complex - 40, ' more'))
    }
    msg = paste0(num_complex, " VCFs had some calls with unsupported ALT alleles(look up the VCF spec for possible explanations of these): ",
                 paste(with_complex, collapse = ', '), '. The calls were not included in MAF output.')
    warning(pretty_message(msg, emit = F))
  }
  maf[, c("possible_multi_sample", "has_complex") := NULL]
  
  
  no_calls = setdiff(sample_ids, maf$Tumor_Sample_Barcode)
  if(length(no_calls) > 0) {
    pretty_message(paste0("FYI, some samples had no variant calls: ", paste(no_calls, collapse = ', '), '.'))
  }
  
  setcolorder(maf, c("Chromosome", "Start_Position", "Reference_Allele", "Tumor_Allele", "Tumor_Sample_Barcode", "FILTER"))
  return(maf[])
}

#' Internal VCF parser
#' 
#' Used by vcfs_to_maf_table()
#' 
#' @param vcf VCF filename or VCF-like data.table.
#' @param sample_id 1-length sample identifier.
#' @param vcf_name 1-length identifier used in some user messages.
#' @return MAF-like data.table
#' @keywords internal
read_vcf = function(vcf, sample_id, vcf_name = sample_id) {
  if(is.character(vcf) && length(vcf) == 1) {
    vcf = tryCatch(fread(vcf, skip = '#CHROM'),
                   error = function(e) {
                     if(conditionMessage(e) %like% "skip='#CHROM' not found")  {
                       msg = paste0("Could not find header line for ", vcf_name, ". It should start with '#CHROM'.")
                       stop(pretty_message(msg, emit = F))
                     } else {
                       stop(e)
                     }
                   })
  } else if(is(vcf, "data.frame")) {
    vcf = as.data.table(vcf)
    # match the leading hash that would appear in a VCF file
    if('CHROM' %in% names(vcf) & ! '#CHROM' %in% names(vcf)) {
      setnames(vcf, 'CHROM', '#CHROM')
    }
  } else {
    stop("vcf should be the path to a VCF file or a data.table with VCF-style columns.")
  }
  
  if(! is.character(sample_id) || length(sample_id) != 1) {
    stop('sample_id should be 1-length character (this function only supports single-tumor-sample VCFs).')
  }
  required_vcf_cols =  c('#CHROM', 'POS', 'REF', 'ALT')
  other_vcf_cols = c('ID', 'QUAL', 'FILTER', 'INFO', 'FORMAT') # currently only checking FILTER
  missing_required = setdiff(required_vcf_cols, names(vcf))
  if(length(missing_required) > 0) {
    stop('Required VCF columns are missing from ', vcf_name, ': ', paste(missing_required, collapse = ', '), '.')
  }
  
  if(length(names(vcf)[names(vcf) %in% c(required_vcf_cols, 'FILTER')]) != length(required_vcf_cols) + 
     ifelse('FILTER' %in% names(vcf), 1, 0)) {
    stop("Column names are repeated in ", vcf_name, '.')
  }
  
  if(vcf[, .N] == 0) {
    return(data.table())
  }
  missing_other = setdiff(other_vcf_cols, names(vcf))
  if('FILTER' %in% missing_other) {
    vcf[, FILTER := NA_character_]
  }

  vcf_sample_cols = setdiff(names(vcf), c(required_vcf_cols, other_vcf_cols))
  
  vcf = vcf[, .(Chromosome = as.character(`#CHROM`), Start_Position = as.numeric(POS), 
                Reference_Allele = toupper(REF), Tumor_Allele = toupper(ALT), FILTER)]
  vcf[, possible_multi_sample := length(vcf_sample_cols) > 2]
  
  if(anyNA(vcf[, .(Chromosome, Start_Position, Reference_Allele, Tumor_Allele)])) {
    stop("There are NA values in CHROM/POS/REF/ALT of ", vcf_name, '.')
  }
  
  # split records with multiple ALT alleles into separate rows
  which_multi_allele = vcf[Tumor_Allele %like% ',', which = T]
  if(length(which_multi_allele) > 0) {
    vcf[, rn := 1:.N]
    to_split = vcf[which_multi_allele]
    alt_by_rn = to_split[, strsplit(Tumor_Allele, ','), by = 'rn'][, .(rn, Tumor_Allele = V1)]
    to_split[, Tumor_Allele := NULL]
    splitted = to_split[alt_by_rn, on = 'rn']
    
    # After splitting alleles, some reference/tumor alleles are longer than they need to be.
    # For example, TAA->T,TA gets split to TAA->T and TAA->TA. The latter should be AA->A (with POS incremented).
    # Iteratively remove the first 5' base from REF/ALT as long as the first two bases match.
    # (This preserves the one shared base that is kept in VCF format.)
    while(TRUE) {
      splitted[,  needs_left_trim := mapply(function(ref, alt) { nchar(alt) > 1 & ref %like% paste0('^', substr(alt, 1, 2))}, 
                                            Reference_Allele, Tumor_Allele) ]
      if(! any(splitted$needs_left_trim)) {
        splitted[, needs_left_trim := NULL]
        break
      }
      splitted[needs_left_trim == TRUE, c("Start_Position", "Reference_Allele", "Tumor_Allele") := 
                 .(Start_Position + 1, substr(Reference_Allele, 2, nchar(Reference_Allele)),
                   substr(Tumor_Allele, 2, nchar(Tumor_Allele)))]
    }
    
    # Matching 3' bases should be dropped since they are logically the same variants
    # (exception forthcoming). Examples: AAT>CAT should be A>C; AAT>AGCAAT should be
    # AA>AGCAA; AAAC>AC should be AAA>A. Exception: We do not include the leftmost base
    # (which is also 3' when length = 1) due to inclusion of flanking reference base. That
    # is, AA>A and A > AA should not be trimmed. Therefore, we drop 3' bases iteratively.
    while(TRUE) {
      splitted[, needs_right_trim := mapply(
        function(ref, alt) {
          # Avoid using any 5' base in the comparison by requiring nchar > 1.
          nchar(ref) > 1 & nchar(alt) > 1 &substr(ref, nchar(ref), nchar(ref)) == substr(alt, nchar(alt), nchar(alt))
        }, Reference_Allele, Tumor_Allele)]
      if(! any(splitted$needs_right_trim)) {
        splitted[, needs_right_trim := NULL]
        break
      }
      splitted[needs_right_trim == TRUE, c("Reference_Allele", "Tumor_Allele") := 
                 .(substr(Reference_Allele, 1, nchar(Reference_Allele) -1),
                   substr(Tumor_Allele, 1, nchar(Tumor_Allele) - 1))]
    }
    
    vcf = vcf[!which_multi_allele]
    vcf = rbind(vcf, splitted[, .SD, .SDcols = names(vcf)], fill = T)[order(rn), -"rn"]
  }
  
  which_complex = vcf[! Reference_Allele %like% '^[ACTG]+$' | ! Tumor_Allele %like% '^[ACTG]+$', which = TRUE]
  num_complex = length(which_complex)
  vcf[, has_complex := num_complex > 0]
  if(num_complex > 0) {
    vcf = vcf[! which_complex]
  }
  
  # To avoid confusion, remove non-variant records.
  which_not_variant = vcf[Reference_Allele == Tumor_Allele, which = T]
  num_not_variant = length(which_not_variant)
  if(num_not_variant > 0) {
    warning(num_not_variant, ' VCF records in ', vcf_name, ' have REF=ALT. These will not be included in output.')
    vcf = vcf[! which_not_variant]
  }
  
  # VCF insertions list the preceding reference base and repeat it in the ALT allele, like T>TGG.
  # Replace the REF allele with "-" and strip the leading reference base from the ALT allele.
  # Note that Start_Position doesn't change.
  vcf[nchar(Reference_Allele) == 1 & substr(Reference_Allele, 1, 1) == substr(Tumor_Allele, 1, 1),
      c("Reference_Allele", "Tumor_Allele") := .('-', substr(Tumor_Allele, 2, nchar(Tumor_Allele)))]
  
  # Deletions, like "GGGA">"G", include the preceding reference base in both the REF and ALT fields.
  # Remove the leading reference base, increment Start_Position, and change Tumor_Allele to '-' to get something like "GGA">"-"
  # Only deletions should have REF[1] and ALT[1] matching.
  vcf[substr(Reference_Allele, 1, 1) == substr(Tumor_Allele, 1, 1), 
      c("Start_Position", "Reference_Allele", "Tumor_Allele") := .(Start_Position + 1, 
                                                                   substr(Reference_Allele, 2, nchar(Reference_Allele)),
                                                                   '-')]
  
  vcf[, Tumor_Sample_Barcode := sample_id]
  return(vcf[])
}
