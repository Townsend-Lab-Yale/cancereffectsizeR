#' Find likely MNVs in an MAF table
#' 
#' Same-sample variants with 2bp of other variants are likely MNVs.
#' 
#' @param maf a valid MAF-style data.table
#' @return a table with MAF records 
#' @keywords internal
detect_mnv = function(maf) {
  # MNVs are only possible in sample/chromosome combinations with more than one MAF record
  poss_mnv = maf[order(Unique_Patient_Identifier, Chromosome, Start_Position)][, .SD[.N > 1], 
                                                                               by = c("Unique_Patient_Identifier", "Chromosome")]
  
  if(poss_mnv[, .N] == 0) {
    return(poss_mnv)
  }
  poss_mnv[, dist_to_prev := c(Inf, diff(Start_Position)), by = c("Unique_Patient_Identifier", "Chromosome")]
  poss_mnv[, dist_to_next := c(dist_to_prev[2:.N], Inf), by = c("Unique_Patient_Identifier", "Chromosome")]
  poss_mnv[dist_to_prev < 3 | dist_to_next < 3, is_mnv := T]
  
  # organize into groups of likely multi-nucleotide events
  mnv = poss_mnv[is_mnv == T]
  mnv[, start_of_group := dist_to_prev > 2]
  mnv[, mnv_group := cumsum(start_of_group)]
  return(mnv)
}
