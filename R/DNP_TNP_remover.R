#' DNP remover
#'
#' Function that removes likely dinucleotide variants that are actually labeled as single nucleotide variants.
#'  Also removes nucleotides 2 positions apart, as these are likely not single nucleotide events
#'  (analysis reveals that read counts for nucleotides two positions apart are almost
#'  perfectly correlated, meaning it is likely sequencing error).
#'  If the tumor sample names have TCGA naming convention and are
#'  under a column header named Tumor_Sample_Barcode,
#'  you can also automatically remove recurrent tumors (and get a warning
#'  about tumor types other than primary)
#'
#' @param MAF MAF file, with relevant fields consisten with the GDC MAF specification (\url{https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/}).
#' @param delete_recur T/F statement about whether to remove recurrent (02) tumors.
#' @export



# Function that removes DNP that are actually labeled as single nucleotide variants.
# Also removes nucleotides 2 positions apart, as these are likely not single nucleotide events
# (analysis reveals that read counts for nucleotides two positions apart are almost
# perfectly correlated, meaning it is likely sequencing error)


DNP_TNP_remover <- function(MAF,delete_recur=F){
  message("Removing possible DNP")
  remove_it <- MAF
  final_MAF <- NULL
  tumor.list <- unique(remove_it$Unique_patient_identifier)
  counter <- 0
  possible.DNP <- NULL
  for(i in 1:length(tumor.list)){
    to.delete <- NULL
    this.tumor_full <- remove_it[which(remove_it$Unique_patient_identifier==tumor.list[i] & remove_it$Variant_Type!="SNP"),]
    this.tumor <- remove_it[which(remove_it$Unique_patient_identifier==tumor.list[i] & remove_it$Variant_Type=="SNP"),]
    for(j in 1:nrow(this.tumor)){
      if(length(which(this.tumor$Chromosome==this.tumor$Chromosome[j] &
                      (this.tumor$Start_Position==this.tumor$Start_Position[j]+1 |
                       this.tumor$Start_Position==this.tumor$Start_Position[j]-1 )))>0){

        to.delete <- c(to.delete,j)
        counter <- counter+1
      }
    }
    for(j in 1:nrow(this.tumor)){
      if(length(which(this.tumor$Chromosome==this.tumor$Chromosome[j] &
                      (this.tumor$Start_Position==this.tumor$Start_Position[j]+2 |
                       this.tumor$Start_Position==this.tumor$Start_Position[j]-2 )))>0){

        to.delete <- c(to.delete,j)
        counter <- counter+1
      }
    }
    if(length(to.delete)>0){
      this.tumor_full <- rbind(this.tumor[-to.delete,],this.tumor_full)
      final_MAF <- rbind(final_MAF,this.tumor_full)
    }else{
      this.tumor_full <- rbind(this.tumor,this.tumor_full)
      final_MAF <- rbind(final_MAF,this.tumor_full)
    }
  }
  message(paste("Total count of potential DNP removed: ", counter))
  message("DNP removal complete")

  # If we only want primary for this analysis
  if(delete_recur){
    message("Deleting any mutations detected in TCGA recurrent tumors")

    if("Tumor_Sample_Barcode" %in% colnames(final_MAF)){
      tumors.unique <- unique(final_MAF$Tumor_Sample_Barcode) # get all the tumor names
      if(length(which(startsWith(x = tumors.unique,prefix = "TCGA")))>1){

        tumors.unique <- tumors.unique[startsWith(x = tumors.unique,prefix = "TCGA")] # get all the tumor names that start with "TCGA"

        tumors.split <- strsplit(tumors.unique,split = "")
        holder.vec <- rep(NA,length(tumors.split))
        for(i in 1:length(tumors.split)){
          holder.vec[i] <- paste(tumors.split[[i]][14:15],sep="",collapse = "")
        }

        if(any(unique(holder.vec)!="01")){ #If any do not equal 01, primary tumor
          warning(paste("Tumor sample barcode sample types exist in these TCGA data besides 01! Removing all 02 "))


          if(length(which(holder.vec=="02"))>0){
            final_MAF <- final_MAF[-which(final_MAF$Tumor_Sample_Barcode %in%  tumors.unique[which(holder.vec=="02")]),]
          }
        }



      }
    }
  }

  return(final_MAF)
}
