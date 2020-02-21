
## Build codon_point_mutation_dict
## for each codon, find all codons that can be generated via one point mutation to the original codon
## keys = original codon; values = list where each name is a new codon, value is the resulting amino acid 

# get all 64 possible codons
all_codons = character(64)
bases = c("A", "T", "C", "G")
j = 1
for (nt1 in bases) {
  for(nt2 in bases) {
    for(nt3 in bases) {
      all_codons[j] = paste0(nt1, nt2, nt3)
      j = j+1
    }
  }
}


# key = original codon (as string); value = list where each name is a new codon, value is the resulting amino acid 
codon_point_mutation_dict = new.env()

for (codon in all_codons) {
  codon = DNAString(codon)
  poss_mutated_codons = list()
  # replace each base with all 4 bases (will remove the non-mutated codons at the end)
  j = 1
  for (nt in bases) {
    for (pos in 1:3) {
      poss_mutated_codons[[j]] =  Biostrings::replaceAt(codon, IRanges::IRanges(pos, pos), nt)
      j = j + 1
    }
  }
  # remove the 3 unmutated codons that were produced
  poss_mutated_codons = DNAStringSet(poss_mutated_codons)
  poss_mutated_codons = poss_mutated_codons[poss_mutated_codons != codon]
  poss_mutated_codons = sapply(poss_mutated_codons, as.character) # done with DNAStrings
  poss_amino_acids = sapply(poss_mutated_codons, function(x) as.character(AA_translations[AA_translations$Nucs == x, "AA_short"]))
  codon_point_mutation_dict[[as.character(codon)]] = setNames(as.list(poss_amino_acids), poss_mutated_codons)
}
usethis::use_data(codon_point_mutation_dict, internal = TRUE)
