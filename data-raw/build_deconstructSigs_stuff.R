# The 96 trinuc-context-specific SNVs, in the order used by deconstructSigs data structures
# Note that the central reference trinucleotides are always C/T, as A/G reference are 
# represented using the reverse complement trinucleotide context
build_deconstructSigs_trinuc_string = function() {
	return(c("A[C>A]A",
	  		 "A[C>A]C",
	  		 "A[C>A]G",
	  		 "A[C>A]T",
	  		 "C[C>A]A",
	  		 "C[C>A]C",
	  		 "C[C>A]G",
	  		 "C[C>A]T",
	  		 "G[C>A]A",
	  		 "G[C>A]C",
	  		 "G[C>A]G",
	  		 "G[C>A]T",
	  		 "T[C>A]A",
	  		 "T[C>A]C",
	  		 "T[C>A]G",
	  		 "T[C>A]T",
	  		 "A[C>G]A",
	  		 "A[C>G]C",
	  		 "A[C>G]G",
	  		 "A[C>G]T",
	  		 "C[C>G]A",
	  		 "C[C>G]C",
	  		 "C[C>G]G",
	  		 "C[C>G]T",
	  		 "G[C>G]A",
	  		 "G[C>G]C",
	  		 "G[C>G]G",
	  		 "G[C>G]T",
	  		 "T[C>G]A",
	  		 "T[C>G]C",
	  		 "T[C>G]G",
	  		 "T[C>G]T",
	  		 "A[C>T]A",
	  		 "A[C>T]C",
	  		 "A[C>T]G",
	  		 "A[C>T]T",
	  		 "C[C>T]A",
	  		 "C[C>T]C",
	  		 "C[C>T]G",
	  		 "C[C>T]T",
	  		 "G[C>T]A",
	  		 "G[C>T]C",
	  		 "G[C>T]G",
	  		 "G[C>T]T",
	  		 "T[C>T]A",
	  		 "T[C>T]C",
	  		 "T[C>T]G",
	  		 "T[C>T]T",
	  		 "A[T>A]A",
	  		 "A[T>A]C",
	  		 "A[T>A]G",
	  		 "A[T>A]T",
	  		 "C[T>A]A",
	  		 "C[T>A]C",
	  		 "C[T>A]G",
	  		 "C[T>A]T",
	  		 "G[T>A]A",
	  		 "G[T>A]C",
	  		 "G[T>A]G",
	  		 "G[T>A]T",
	  		 "T[T>A]A",
	  		 "T[T>A]C",
	  		 "T[T>A]G",
	  		 "T[T>A]T",
	  		 "A[T>C]A",
	  		 "A[T>C]C",
	  		 "A[T>C]G",
	  		 "A[T>C]T",
	  		 "C[T>C]A",
	  		 "C[T>C]C",
	  		 "C[T>C]G",
	  		 "C[T>C]T",
	  		 "G[T>C]A",
	  		 "G[T>C]C",
	  		 "G[T>C]G",
	  		 "G[T>C]T",
	  		 "T[T>C]A",
	  		 "T[T>C]C",
	  		 "T[T>C]G",
	  		 "T[T>C]T",
	  		 "A[T>G]A",
	  		 "A[T>G]C",
	  		 "A[T>G]G",
	  		 "A[T>G]T",
	  		 "C[T>G]A",
	  		 "C[T>G]C",
	  		 "C[T>G]G",
	  		 "C[T>G]T",
	  		 "G[T>G]A",
	  		 "G[T>G]C",
	  		 "G[T>G]G",
	  		 "G[T>G]T",
	  		 "T[T>G]A",
	  		 "T[T>G]C",
	  		 "T[T>G]G",
	  		 "T[T>G]T"))
}

# build the data frame translating trinuc-context SNV mutations into the format used by deconstructSigs
# e.g.: a C>T where in the reference sequence, C is in an ACA context (i.e., ACA:T) becomes "A[C>T]A"
# when reference is A or G, deconstructSigs represents with the reverse complement (e.g., AGT:A -> "A[C>T]T")
build_trinuc_translator = function(deconstructSigs_trinuc_string) {
	# create data frame where each row lists one possible combination of upstream, ref, downstream, and mutation
	nts = c("A", "T", "G", "C")

	all_trinucs <- expand.grid(nts, nts, nts, nts,stringsAsFactors = F)
	colnames(all_trinucs) <- c("upstream","ref","downstream","mut")

	# drop the rows where reference matches variant
	all_trinucs <- all_trinucs[-which(all_trinucs$ref == all_trinucs$mut),]

	# add a column with the "deconstructSigs" mutation format (e.g., A[C>T]G)
	all_trinucs$deconstructSigs_format <- NA
	for(i in 1:nrow(all_trinucs)){
	  if(all_trinucs[i,"ref"] %in% c("C","T")){
	    all_trinucs$deconstructSigs_format[i] <- paste(all_trinucs$upstream[i],"[",all_trinucs$ref[i],">",
	    											   all_trinucs$mut[i],"]",all_trinucs$downstream[i], sep = "")
	  }else{
	  	# reverse compelement of mutations to A and G
	    all_trinucs$deconstructSigs_format[i] <- paste(toupper(seqinr::comp(all_trinucs$downstream[i])),
	                                                   "[",toupper(seqinr::comp(all_trinucs$ref[i])),">",
	                                                   toupper(seqinr::comp(all_trinucs$mut[i])),"]",
	                                                   toupper(seqinr::comp(all_trinucs$upstream[i])),sep="")
	  }
	}
	rownames(all_trinucs) <- paste(all_trinucs$upstream,all_trinucs$ref,all_trinucs$downstream,":",all_trinucs$mut,sep="")
	return(all_trinucs)
}
