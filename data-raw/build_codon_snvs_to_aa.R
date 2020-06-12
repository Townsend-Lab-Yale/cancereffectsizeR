library(Biostrings)

# poss_mut: A reference list that, given a particular genomic coding trinucleotide sequence,
# and an amino acid of interest, gives all possible SNVs that would generate the amino acid.
# For each trinucleotide sequence any amino acid, there is a three-item list that gives all
# substitions in the first, second, and third positions that produce the aminoacid.

# Examples:
# poss_mut$TAG$Trp: TAG (which encodes a stop codon) instead encodes Trp after an A->G substitution at position 2
# poss_mut$GTA$Val: GTA (Valine) continues to code Valine regardless of the the third nucleotide, so
# poss_mut$GTA$Val[[3]] returns c("C", "G", "T"). "A" is not included because it is not a substitution.
# poss_mut$AGC$Arg: Shows that substitutions in first and third positions can cause Ser -> Arg


build_codon_snvs_to_aa = function() {
  # a handy conversion from short AA abbreviation to long (e.g., K -> Lys)
  aa_long = AMINO_ACID_CODE
  aa_long["*"] = "STOP" # add in stop codon
  
  
  poss_mut = list()
  for (nt1 in DNA_BASES) {
    for (nt2 in DNA_BASES) {
      for (nt3 in DNA_BASES) {
        # for the given three nts, construct a codon and get its amino acid
        codon = DNAString(xscat(nt1, nt2, nt3)) # need to cast BString as DNAString for translation
        amino_acid = translate(codon)
        codon_str = as.character(codon)
        
        # make empty list structure to hold possible mutations
        poss_mut[[codon_str]] = list()
        aa_code = aa_long[c(AA_STANDARD, "*")] # get long abbreviations of 20 amino acids, plus stop
        for (i in aa_code) {
          poss_mut[[codon_str]][[i]] = list()
          for (j in 1:3) {
            poss_mut[[codon_str]][[i]][[j]] = character()
          }
        }
        
        # see ?Biostrings::translate, but if you don't use no.init.codon, you'll get some mistranslations (e.g., CTG as M instead of L)
        for (mut in setdiff(DNA_BASES, nt1)) {
          new_aa = as.character(translate(replaceLetterAt(codon, 1, mut), no.init.codon = T))
          new_aa = aa_long[new_aa]
          poss_mut[[codon_str]][[new_aa]][[1]] = c(poss_mut[[codon_str]][[new_aa]][[1]], mut)
        }
        for (mut in setdiff(DNA_BASES, nt2)) {
          new_aa = as.character(translate(replaceLetterAt(codon, 2, mut), no.init.codon = T))
          new_aa = aa_long[new_aa]
          poss_mut[[codon_str]][[new_aa]][[2]] = c(poss_mut[[codon_str]][[new_aa]][[2]], mut)
        }
        for (mut in setdiff(DNA_BASES, nt3)) {
          new_aa = as.character(translate(replaceLetterAt(codon, 3, mut), no.init.codon = T))
          new_aa = aa_long[new_aa]
          poss_mut[[codon_str]][[new_aa]][[3]] = c(poss_mut[[codon_str]][[new_aa]][[3]], mut)
        }
      }
    }
  }
  return(poss_mut)
}

