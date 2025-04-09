build_two_codon_dbs_key = function() {
  # seqinr wants lower case ¯\_(ツ)_/¯ 	
  nt = c('a', 'c', 'g', 't')
  
  two_dna_codon = as.data.table(expand.grid(a = nt, b = nt, c = nt, d = nt, e = nt, f = nt,
                                            stringsAsFactors = FALSE))
  two_dna_codon[, aa1 := apply(.SD, 1, seqinr::translate), .SDcols = 1:3]
  two_dna_codon[, aa2 := apply(.SD, 1, seqinr::translate), .SDcols = 4:6]
  
  # Back to uppercase.
  two_dna_codon[, names(.SD) := lapply(.SD, toupper)]
  two_dna_codon[, double_aa := paste0(aa1, aa2)]
  
  two_dna_codon[, spanning_dinuc := paste0(c, d)]
  two_dna_codon[, rest_of_seq := paste0(a, b, e, f)]
  
  by_rest_of_seq = split(two_dna_codon[, .(double_aa, spanning_dinuc)], two_dna_codon$rest_of_seq)
  two_codon_dbs_key = lapply(by_rest_of_seq, function(x) {
    by_aa = x[, .(dn = list(spanning_dinuc)), by = 'double_aa']
    return(setNames(by_aa$dn, by_aa$double_aa))
  })
  return(two_codon_dbs_key)
}

