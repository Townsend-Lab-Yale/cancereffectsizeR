#'  Remove likely multi-nucleotide mutation events
#'  
#'  Function that removes DNP that are actually labeled as single nucleotide variants.
#   Also removes nucleotides 2 positions apart, as these are likely not single nucleotide events
#   (analysis reveals that read counts for nucleotides two positions apart are almost
#   perfectly correlated, meaning it is likely sequencing error)
#'
#' @param MAF MAF file, with relevant fields consisten with the GDC MAF specification (\url{https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/}).
#' @param delete_recur T/F statement about whether to remove recurrent (02) tumors.
#' @export
#' @keywords internal
DNP_TNP_remover <- function(MAF,delete_recur=F){
  #message("Removing possible DNP and TNP")
  # sort by Tumor, Chromosome, and then position
  #  ordering with this syntax so rows get sorted as they would on a df (to keep a test functional)
  ordering = with(MAF, order(Unique_Patient_Identifier,Chromosome,Start_Position))
  MAF <- MAF[ordering,]
  MAF$difference <- c(NA,diff(MAF$Start_Position))

  MAF$same_tumor_chrom <- c(NA,((MAF[2:nrow(MAF),Chromosome]==MAF[1:(nrow(MAF)-1),Chromosome]) &
                                  (MAF[2:nrow(MAF),Unique_Patient_Identifier]==MAF[1:(nrow(MAF)-1),Unique_Patient_Identifier])))

  DNP_and_TNP <- which(MAF$difference %in% c(1,2) & MAF$same_tumor_chrom==T)
  DNP_and_TNP_prior <- DNP_and_TNP - 1


  to_remove <- unique(c(DNP_and_TNP,DNP_and_TNP_prior))

  dnp_and_tnp = MAF[to_remove, 1:5] # to-do: make this less cryptic
  dnp_and_tnp = dnp_and_tnp[with(dnp_and_tnp, order(Chromosome, Start_Position, Unique_Patient_Identifier)),]

  if(length(to_remove) != 0) {
    MAF <- MAF[-to_remove,]
  }
  
  


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

  return(list(kept = MAF, removed = dnp_and_tnp))
}
