#' Prep ASCAT3 copy segment data
#' 
#' @param segments docs coming
#' @param refset reference data package (or custom reference data)
#' @export
prep_ASCAT3_segments = function(segments, refset) {
  if(is.data.frame(segments)) {
    segments = as.data.table(segments)
  } else if(rlang::is_scalar_character(segments)) {
    if(file.exists(segments)) {
      segments = fread(segments)
    } else {
      stop('File ', segments, ' not found.')
    }
  } else if(! is.data.table(segments)) {
    stop('segments should be a file path or data.table.')
  }
  
  # To-do: validate ASCAT3 calls
  cna_calls = copy(segments)
  req_cols = c('patient_id', 'chr', 'start', 'end', 'total_copy', 'nMajor', 'nMinor', 'sex')
  
  if(! all(req_cols %in% names(cna_calls))) {
    # Note two columns stay the same
    setnames(cna_calls, 
             c('patient_id', 'Chromosome', 'Start', 'End', 'Copy_Number', 'Major_Copy_Number', 'Minor_Copy_Number', 'sex'),
             req_cols)
  }
  
  cna_calls = cna_calls[, .SD, .SDcols = req_cols]
  
  if(anyNA(cna_calls$sex)) {
    stop('NA values in sex column not yet supported.')
  }

  if(is.environment(refset)) {
    refset = as.character(substitute(refset))
  }
  refset_info = parse_refset_name(refset)
  #coordinates = get_ref_data(refset_info$)
  
  coordinates = get_ref_data(refset_info$data_dir, 'arm_coordinates')
  genome_info = get_ref_data(refset_info$data_dir, 'genome_build_info')
  get_ref_data(refset_info$data_dir, 'arm_coordinates')
  setnames(cna_calls, 'patient_id', 'sample')
  prepped_calls = .prep_ascat_segments(calls = cna_calls, chr_coordinates = coordinates, 
                                      genome_build = genome_info$build_name)
  return(prepped_calls)
}


