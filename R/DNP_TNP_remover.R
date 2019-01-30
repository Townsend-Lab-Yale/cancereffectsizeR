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
  message("Removing possible DNP and TNP")
  # sort by Tumor, Chromosome, and then position
  MAF <- MAF[with(MAF, order(Unique_patient_identifier,Chromosome,Start_Position)),]

  MAF$difference <- c(NA,diff(MAF[,"Start_Position"]))

  MAF$same_tumor_chrom <- c(NA,((MAF[2:nrow(MAF),"Chromosome"]==MAF[1:(nrow(MAF)-1),"Chromosome"]) &
                                  (MAF[2:nrow(MAF),"Unique_patient_identifier"]==MAF[1:(nrow(MAF)-1),"Unique_patient_identifier"])))

  DNP_and_TNP <- which(MAF$difference %in% c(1,2) & MAF$same_tumor_chrom==T)
  DNP_and_TNP_prior <- DNP_and_TNP - 1


  to_remove <- unique(c(DNP_and_TNP,DNP_and_TNP_prior))

  if(length(to_remove)>0){
    MAF <- MAF[-to_remove,]
  }
  # remove_it <- MAF
  # final_MAF <- NULL
  # tumor_list <- unique(remove_it$Unique_patient_identifier)
  # counter <- 0
  # possible_DNP <- NULL
  # for(i in 1:length(tumor_list)){
  #   to_delete <- NULL
  #   this_tumor_full <- remove_it[which(remove_it$Unique_patient_identifier==tumor_list[i] & remove_it$Variant_Type!="SNP"),]
  #   this_tumor <- remove_it[which(remove_it$Unique_patient_identifier==tumor_list[i] & remove_it$Variant_Type=="SNP"),]
  #   for(j in 1:nrow(this_tumor)){
  #     if(any(this_tumor$Chromosome==this_tumor$Chromosome[j] &
  #                     (this_tumor$Start_Position==this_tumor$Start_Position[j]+1 |
  #                      this_tumor$Start_Position==this_tumor$Start_Position[j]-1 |
  #                      this_tumor$Start_Position==this_tumor$Start_Position[j]+2 |
  #                      this_tumor$Start_Position==this_tumor$Start_Position[j]-2 ), na.rm=TRUE)){
  #
  #       to_delete <- c(to_delete,j)
  #       counter <- counter+1
  #     }
  #   }
  #
  #   if(length(to_delete)>0){
  #     this_tumor_full <- rbind(this_tumor[-to_delete,],this_tumor_full)
  #     final_MAF <- rbind(final_MAF,this_tumor_full)
  #   }else{
  #     this_tumor_full <- rbind(this_tumor,this_tumor_full)
  #     final_MAF <- rbind(final_MAF,this_tumor_full)
  #   }
  # }
  message(paste("Total count of potential DNP removed: ", length(to_remove)))
  message("DNP and TNP removal complete")

  # If we only want primary for this analysis
  if(delete_recur){
    message("Deleting any mutations detected in TCGA recurrent tumors")

    if("Tumor_Sample_Barcode" %in% colnames(MAF)){
      tumors_unique <- unique(MAF$Tumor_Sample_Barcode) # get all the tumor names
      if(sum(startsWith(x = tumors_unique,prefix = "TCGA"))>1){

        tumors_unique <- tumors_unique[startsWith(x = tumors_unique,prefix = "TCGA")] # get all the tumor names that start with "TCGA"

        holder_vec <- substr(tumors_unique, 14, 15)

        if(any(unique(holder_vec)!="01")){ #If any do not equal 01, primary tumor
          warning(paste("Tumor sample barcode sample types exist in these TCGA data besides 01! Removing all 02 "))

          if(any(holder_vec=="02", na.rm=TRUE)){
            MAF <- MAF[!(MAF$Tumor_Sample_Barcode %in% tumors_unique[which(holder_vec=="02")]),]
          }
        }
      }
    }
  }

  return(MAF)
}
