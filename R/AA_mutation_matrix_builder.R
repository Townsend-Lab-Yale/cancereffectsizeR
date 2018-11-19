# making a list of all possible trinucleotides to call for each amino acid variant.

library(tidyverse)
library(cancereffectsizeR)

load("data/deconstructSigs_trinuc_string.RData")

all.nucs <- expand.grid(c("A","T","G","C"),c("A","T","G","C"),c("A","T","G","C"),c("A","T","G","C"),c("A","T","G","C"),stringsAsFactors = F)

all.nucs$combination <- apply(all.nucs,1,function(x)paste(x,collapse=""))

AA_mutation_list <- vector("list",nrow(all.nucs))
names(AA_mutation_list) <- all.nucs$combination

AA_translations$AA_short <- as.character(AA_translations$AA_short)
AA_translations$Nucs <- as.character(AA_translations$Nucs)
AA_translations$AA_name <- as.character(AA_translations$AA_name)
AA_translations$AA_letter <- as.character(AA_translations$AA_letter)
rownames(AA_translations) <- AA_translations$Nucs

nucs <- c("A","T","G","C")


for(i in 1:nrow(all.nucs)){
  this_codon <- as.character(all.nucs[i,c("Var2","Var3","Var4")])
  this_5seq <- as.character(all.nucs[i,c("Var1","Var2","Var3","Var4","Var5")])
  this_amino_acid <- AA_translations[paste(this_codon,collapse = ""),"AA_short"]
  alt_AA=list(Phe= NULL,
              Leu= NULL,
              Ser= NULL,
              Tyr= NULL,
              Cys= NULL,
              Trp = NULL,
              Pro=NULL,
              His=NULL,Gln=NULL,Arg=NULL,Ile=NULL,Met=NULL,Thr=NULL,Asn=NULL,Lys=NULL,Val=NULL,Ala=NULL,Asp=NULL,Glu=NULL,Gly=NULL,STOP=NULL)

  for(codon_pos in 2:4){

    for(alt_allele in 1:3){
      new_codon <- this_codon
      new_codon[codon_pos-1] <- setdiff(y = this_codon[codon_pos-1],x = nucs)[alt_allele]

      new_5seq <- this_5seq
      new_5seq[codon_pos] <- setdiff(y = this_codon[codon_pos-1],x = nucs)[alt_allele]

      if(this_codon[codon_pos-1] %in% c("C","T")){

        this_trinuc <- new_5seq[c(codon_pos-1,codon_pos,codon_pos+1)]

        trinuc_to_paste <- paste(this_trinuc[1],"[",this_codon[codon_pos-1],">",this_trinuc[2],"]",this_trinuc[3],sep = "")

      }

      if(this_codon[codon_pos-1] %in% c("A","G")){

        this_trinuc <- new_5seq[c(codon_pos-1,codon_pos,codon_pos+1)]
        this_trinuc <- rev(toupper(seqinr::comp(this_trinuc)))

        trinuc_to_paste <- paste(this_trinuc[1],"[",toupper(seqinr::comp(this_codon[codon_pos-1])),">",this_trinuc[2],"]",this_trinuc[3],sep = "")

      }

      alt_AA[[AA_translations[paste(new_codon,collapse = ""),"AA_short"]]] <-
        c(alt_AA[[AA_translations[paste(new_codon,collapse = ""),"AA_short"]]],trinuc_to_paste)

    }
  }

  AA_mutation_list[[i]] <- alt_AA
}


save(AA_mutation_list,file = "data/AA_mutation_list.RData")




# AA_mutation_list[[1]] <- list(K=c("AAA","VVV","AAA","XXX"))
#
# AA_mutation_list[["AAAAA"]]["K"]
# typeof(AA_mutation_list[["AAAAA"]]["K"])
# test.matrix <- matrix(nrow=5,ncol=1,data=1:5)
# rownames(test.matrix) <- c("AAA","VVV","XXX","GGG","TTT")
#
# test.matrix[as.character(unlist(AA_mutation_list[["AAAAA"]]["K"])),]
#
#
# test.df <- data.frame(AA=as.character(unique(AA_translations$AA_letter)),muts=NA)
#
# test.df[1,2] <- list(c("AA","TT"))
#
#
#
#
#
#
#
