# create list of all the trinuc numbers in every gene
# library(cancereffectsizeR)
library(dndscv)
data("RefCDS_TP53splice", package="cancereffectsizeR")

gene_list <- sapply(RefCDS, function(x) x$gene_name)

load("data/deconstructSigs_trinuc_string.RData")

gene_trinuc_comp <- vector("list",length(RefCDS))

for(i in 1:length(gene_trinuc_comp)){

  gene_trinuc <- data.frame(count=rep(0,length(deconstructSigs_trinuc_string)))
  rownames(gene_trinuc) <- deconstructSigs_trinuc_string

  gene_trinuc$downstream <- substr(deconstructSigs_trinuc_string,start = 1,stop = 1)
  gene_trinuc$upstream <- substr(deconstructSigs_trinuc_string,start = 7,stop = 7)
  gene_trinuc$ref <- substr(deconstructSigs_trinuc_string,start = 3,stop = 3)
  gene_trinuc$alt <- substr(deconstructSigs_trinuc_string,start = 5,stop = 5)


  for(j in 1:RefCDS[[i]]$CDS_length){

    if(!(strsplit(as.character(RefCDS[[i]]$seq_cds),split = "")[[1]][j] %in% c("C","T"))){
      this_trinuc <- c(
        strsplit(as.character(RefCDS[[i]]$seq_cds1up),split = "")[[1]][j],
        strsplit(as.character(RefCDS[[i]]$seq_cds),split = "")[[1]][j],
        strsplit(as.character(RefCDS[[i]]$seq_cds1down),split = "")[[1]][j]
      )

      this_trinuc <- rev(toupper(seqinr::comp(this_trinuc)))

    }else{
      this_trinuc <- c(
        strsplit(as.character(RefCDS[[i]]$seq_cds1up),split = "")[[1]][j],
        strsplit(as.character(RefCDS[[i]]$seq_cds),split = "")[[1]][j],
        strsplit(as.character(RefCDS[[i]]$seq_cds1down),split = "")[[1]][j]
      )

    }

    gene_trinuc$count[which(gene_trinuc$upstream== this_trinuc[1] &
                              gene_trinuc$ref== this_trinuc[2] &
                              gene_trinuc$downstream == this_trinuc[3])] <-
      gene_trinuc$count[which(gene_trinuc$upstream== this_trinuc[1] &
                                gene_trinuc$ref== this_trinuc[2] &
                                gene_trinuc$downstream == this_trinuc[3])] + 1


  }


  if(RefCDS[[i]]$strand=="-1"){
    seq_cds2up <- paste(toupper(seqinr::comp(strsplit(as.character(rev(unlist(BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens,paste("chr",RefCDS[[i]]$chr,sep=""),start=RefCDS[[i]]$intervals_cds[,1]+2,end = RefCDS[[i]]$intervals_cds[,2]+2,strand="+")))),split="")[[1]])),collapse = "")

    seq_cds2down <- paste(toupper(seqinr::comp(strsplit(as.character(rev(unlist(BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens,paste("chr",RefCDS[[i]]$chr,sep=""),start=RefCDS[[i]]$intervals_cds[,1]-2,end = RefCDS[[i]]$intervals_cds[,2]-2,strand="+")))),split="")[[1]])),collapse = "")

    seq_cds3up <- paste(toupper(seqinr::comp(strsplit(as.character(rev(unlist(BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens,paste("chr",RefCDS[[i]]$chr,sep=""),start=RefCDS[[i]]$intervals_cds[,1]+3,end = RefCDS[[i]]$intervals_cds[,2]+3,strand="+")))),split="")[[1]])),collapse = "")

    seq_cds3down <- paste(toupper(seqinr::comp(strsplit(as.character(rev(unlist(BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens,paste("chr",RefCDS[[i]]$chr,sep=""),start=RefCDS[[i]]$intervals_cds[,1]-3,end = RefCDS[[i]]$intervals_cds[,2]-3,strand="+")))),split="")[[1]])),collapse = "")

  }


  if(RefCDS[[i]]$strand=="1"){
    seq_cds2up <- paste(strsplit(as.character(unlist(BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens,paste("chr",RefCDS[[i]]$chr,sep=""),start=RefCDS[[i]]$intervals_cds[,1]-2,end = RefCDS[[i]]$intervals_cds[,2]-2,strand="+"))),split="")[[1]],collapse = "")

    seq_cds2down <- paste(strsplit(as.character(unlist(BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens,paste("chr",RefCDS[[i]]$chr,sep=""),start=RefCDS[[i]]$intervals_cds[,1]+2,end = RefCDS[[i]]$intervals_cds[,2]+2,strand="+"))),split="")[[1]],collapse = "")

    seq_cds3up <- paste(strsplit(as.character(unlist(BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens,paste("chr",RefCDS[[i]]$chr,sep=""),start=RefCDS[[i]]$intervals_cds[,1]-3,end = RefCDS[[i]]$intervals_cds[,2]-3,strand="+"))),split="")[[1]],collapse = "")

    seq_cds3down <- paste(strsplit(as.character(unlist(BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens,paste("chr",RefCDS[[i]]$chr,sep=""),start=RefCDS[[i]]$intervals_cds[,1]+3,end = RefCDS[[i]]$intervals_cds[,2]+3,strand="+"))),split="")[[1]],collapse = "")


  }


  gene_trinuc_comp[[i]] <- list(gene_name = gene_list[i],gene_trinuc=gene_trinuc,
                                seq_cds2up=seq_cds2up, seq_cds2down=seq_cds2down,
                                seq_cds3up=seq_cds3up, seq_cds3down=seq_cds3down)


  print(i)

}


save(gene_trinuc_comp,file = "data/gene_trinuc_comp.RData")




