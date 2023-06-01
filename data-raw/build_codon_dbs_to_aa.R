build_codon_dbs_to_aa = function() {
  # a handy conversion from short AA abbreviation to long (e.g., K -> Lys)
  aa_long = AMINO_ACID_CODE
  aa_long["*"] = "STOP" # add in stop codon
  
  nt = c('A', 'C', 'T', 'G')
  dinuc = as.data.table(S4Vectors::expand.grid(nt, nt))[, paste0(Var1, Var2)]
  
  
  poss_mut = list()
  for (nt1 in nt) {
    for (nt2 in nt) {
      for (nt3 in nt) {
        # for the given three nts, construct a codon and get its amino acid
        codon = DNAString(xscat(nt1, nt2, nt3)) # need to cast BString as DNAString for translation
        amino_acid = translate(codon)
        codon_str = as.character(codon)
        
        # make empty list structure to hold possible mutations
        poss_mut[[codon_str]] = list()
        aa_code = aa_long[c(AA_STANDARD, "*")] # get long abbreviations of 20 amino acids, plus stop
        for (i in aa_code) {
          poss_mut[[codon_str]][[i]] = list()
          for (j in 1:2) {
            poss_mut[[codon_str]][[i]][[j]] = character()
          }
        }
        
        # see ?Biostrings::translate, but if you don't use no.init.codon, you'll get some mistranslations (e.g., CTG as M instead of L)
        # Specify all dinucleotides that don't count as DBS (due to one or both bases not changing)
        not_dbs = c(paste0(nt1, nt),
                      paste0(nt, nt2))
        original_codon = codon
        for (mut in setdiff(dinuc, not_dbs)) {
          subseq(codon, 1, 2) = DNAString(mut)
          new_aa = as.character(translate(codon, no.init.codon = T))
          new_aa = aa_long[new_aa]
          poss_mut[[codon_str]][[new_aa]][[1]] = c(poss_mut[[codon_str]][[new_aa]][[1]], mut)
        }
        codon = original_codon
        not_dbs = c(paste0(nt2, nt),
                    paste0(nt, nt3))
        for (mut in setdiff(dinuc, not_dbs)) {
          subseq(codon, 2, 3) = DNAString(mut)
          new_aa = as.character(translate(codon, no.init.codon = T))
          new_aa = aa_long[new_aa]
          poss_mut[[codon_str]][[new_aa]][[2]] = c(poss_mut[[codon_str]][[new_aa]][[2]], mut)
        }
      }
    }
  }
  return(poss_mut)
}
