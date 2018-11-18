# making a list of all possible trinucleotides to call for each amino acid variant.

library(tidyverse)
library(cancereffectsizeR)

all.nucs <- expand.grid(c("A","T","G","C"),c("A","T","G","C"),c("A","T","G","C"),c("A","T","G","C"),c("A","T","G","C"),stringsAsFactors = F)

all.nucs$combination <- apply(all.nucs,1,function(x)paste(x,collapse=""))

AA_mutation_list <- vector("list",nrow(all.nucs))
names(AA_mutation_list) <- all.nucs$combination


for(i in 1:nrow(all.nucs)){
  this_codon <- all.nucs[i,c("Var2","Var3","Var4")]
  this_amino_acid <- seqinr::translate(this_codon)
  Phe= NULL
  Leu= NULL
  Ser= NULL
  Tyr= NULL


}




AA_mutation_list[[1]] <- list(K=c("AAA","VVV","AAA","XXX"))

AA_mutation_list[["AAAAA"]]["K"]
typeof(AA_mutation_list[["AAAAA"]]["K"])
test.matrix <- matrix(nrow=5,ncol=1,data=1:5)
rownames(test.matrix) <- c("AAA","VVV","XXX","GGG","TTT")

test.matrix[as.character(unlist(AA_mutation_list[["AAAAA"]]["K"])),]


test.df <- data.frame(AA=as.character(unique(AA_translations$AA_letter)),muts=NA)

test.df[1,2] <- list(c("AA","TT"))







