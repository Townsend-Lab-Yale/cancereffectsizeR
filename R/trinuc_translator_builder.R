library(seqinr)

load("data/deconstructSigs_trinuc_string.RData")
deconstructSigs_trinuc_string

all_trinucs <- expand.grid(c("A","T","G","C"),c("A","T","G","C"),c("A","T","G","C"),c("A","T","G","C"),stringsAsFactors = F)

colnames(all_trinucs) <- c("upstream","ref","downstream","mut")

all_trinucs <- all_trinucs[-which(all_trinucs$ref == all_trinucs$mut),]

all_trinucs$deconstructSigs_format <- NA

for(i in 1:nrow(all_trinucs)){
  if(all_trinucs[i,"ref"] %in% c("C","T")){
    all_trinucs$deconstructSigs_format[i] <- paste(all_trinucs$upstream[i],"[",all_trinucs$ref[i],">",all_trinucs$mut[i],"]",all_trinucs$downstream[i],sep = "")
  }else{
    all_trinucs$deconstructSigs_format[i] <- paste(toupper(seqinr::comp(all_trinucs$downstream[i])),
                                                   "[",toupper(seqinr::comp(all_trinucs$ref[i])),">",toupper(seqinr::comp(all_trinucs$mut[i])),"]",toupper(seqinr::comp(all_trinucs$upstream[i])),sep="")
  }
}


rownames(all_trinucs) <- paste(all_trinucs$upstream,all_trinucs$ref,all_trinucs$downstream,":",all_trinucs$mut,sep="")

trinuc_translator <- all_trinucs

save(trinuc_translator,file = "data/trinuc_translator.RData")
