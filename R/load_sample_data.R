#' Add sample data 
#' 
#' Insert sample-level data into a CESAnalysis samples table.
#' 
#' @param cesa CESAnalysis
#' @param sample_data data.table or data.frame with a Unique_Patient_Identifier column to
#'   match the CESAnalysis samples table, with one row per sample. (It's okay if some
#'   samples aren't present in the table.) For simplicity, columns besides
#'   Unique_Patient_Identifier must not already be present in the CESAnalysis. (Instead,
#'   build a single sample_data table and load it with one call to this function.)
#' @export
load_sample_data = function(cesa, sample_data) {
  if(! is(cesa, "CESAnalysis")) {
    stop("cesa expected to be a CESAnalysis.")
  }
  
  if(cesa@samples[, .N] == 0) {
    stop("There are no samples loaded into the CESAnalysis.")
  }
  if (is(sample_data, "data.frame")) {
    sample_data = as.data.table(sample_data)
  }
  if (! is(sample_data, "data.table")) {
    stop("sample_data must a be a data.table or data.frame.")
  }
  
  if(! 'Unique_Patient_Identifier' %in% names(sample_data)) {
    stop("sample_data must have a Unique_Patient_Identifier column to match the CESAnalysis sample table.")
  }
  
  # Will check here, and then again after subtracting internal columns
  if(length(sample_data) == 1) {
    stop("There are no data columns in the input sample_data.")
  }
  
  if(sample_data[, .N] == 0) {
    stop("sample_data has 0 rows.")
  }
  
  
  # Check for non-internal columns already present from previous load_sample_data() calls.
  # internal_columns similar to sample_table_template but with some past/future columns
  internal_columns = c('coverage', 'covered_regions', 'group', 'gene_rate_grp', 'regional_rate_grp', 'sig_analysis_grp')
  already_in_maf = setdiff(names(cesa@samples), c("Unique_Patient_Identifier", internal_columns))
  
  dup_columns = unique(names(sample_data)[duplicated(names(sample_data))])
  if(length(dup_columns) > 0) {
    stop("Some column names in sample_data are used more than once: ", 
         paste(dup_columns, collapse = ", "))
  }
  
  reused = intersect(names(sample_data), already_in_maf)
  if (length(reused) > 0) {
    msg = paste0("Ignoring sample_data columns already in samples table (", 
                 paste(reused, collapse = ", ") , ").")
    pretty_message(msg, black = F)
    pretty_message("(Run clear_sample_data() first if you want to overwrite these columns.)", black = F)
    sample_data = sample_data[, .SD, .SDcols = setdiff(names(sample_data), reused)]
  }
  
  missing_samples = setdiff(sample_data$Unique_Patient_Identifier, cesa@samples$Unique_Patient_Identifier)
  num_missing = length(missing_samples)
  if (num_missing > 0) {
    if (num_missing < 25) {
      stop(num_missing, " samples in sample_data are not present in the CESAnalysis sample table: ",
           paste(missing_samples, collapse = ", "))
    } else {
      missing_samples = missing_samples[1:10]
      stop(num_missing, " samples in sample_data are not present in the CESAnalysis sample table, including ",
           paste(missing_samples[1:10], collapse = ", "))
    }
  }
  
  repeated_samples = unique(sample_data$Unique_Patient_Identifier[duplicated(sample_data$Unique_Patient_Identifier)])
  if(length(repeated_samples) > 0) {
    stop("Some samples appear more than once in sample_data: ", 
         paste(repeated_samples, collapse = ", "))
  }
  
  
  # Verify sample_data doesn't contain columns reserved by cancereffectsizeR.
  # Exception: We'll allow the reserved columns if their values exactly match sample table.
  # Also saving regional_rate_group for future use.
  disallowed_columns = intersect(internal_columns, names(sample_data))
  if(length(disallowed_columns) > 0) {
    if (all(disallowed_columns %in% names(cesa@samples)) && 
            identical(all.equal(cesa@samples[sample_data$Unique_Patient_Identifier, ..disallowed_columns, on = "Unique_Patient_Identifier"], 
                  sample_data[, ..disallowed_columns], check.attributes = FALSE), TRUE)) {
        sample_data = sample_data[, ! ..disallowed_columns]
        
        # If just one column, must be Unique_Patient_Identifier with no data columns
        if(length(sample_data) == 1) {
          stop("There are no new data columns in sample_data.")
        }
    } else {
      stop("sample_data contains columns reserved for internal use, and the data don't match. Please remove or rename them:\n",
           paste(disallowed_columns, collapse = ", "))
    }
  }
  
  
  
  
  
  cesa = copy_cesa(cesa)
  cesa = update_cesa_history(cesa, match.call())
  
  # Automatically sorts/keys on Unique_Patient_Identifier, too
  cesa@samples = merge.data.table(cesa@samples, sample_data, by = 'Unique_Patient_Identifier', all.x = T)
  return(cesa)
}


#' Clear sample data
#' 
#' Remove data columns by name from CESAnalysis sample table. You can't clear
#' cancereffectsizeR-generated columns, such as coverage.
#' 
#' @param cesa CESAnalysis
#' @param cols names of data columns to clear
#' @export
#' 
clear_sample_data = function(cesa, cols) {
  if(! is(cesa, "CESAnalysis")) {
    stop("Expected cesa to be CESAnalysis.")
  }
  if(cesa@samples[, .N] == 0) {
    stop("The CESAnalysis samples table is empty.")
  }
  
  if(! is.character(cols) || length(cols) == 0) {
    stop("cols should be a character vector of column names.")
  }
  
  cols = unique(cols)
  
  internal_columns = c('coverage', 'covered_regions', 'group', 'gene_rate_grp', 'regional_rate_grp', 'Unique_Patient_Identifier')
  protected_cols = intersect(cols, internal_columns)
  if (length(protected_cols) > 0) {
    stop("Internal columns can't be cleared: ", paste(protected_cols, collapse = ", "), ".")
  }
  
  missing_cols = setdiff(cols, names(cesa@samples))
  if(length(missing_cols) > 0) {
    stop("Not all specified columns are present in the sample table: ", paste(missing_cols, collapse = ", "), ".")
  }
  
  cesa = copy_cesa(cesa)
  cesa@samples[, (cols) := NULL]
  cesa = update_cesa_history(cesa, match.call())
  return(cesa)
}
