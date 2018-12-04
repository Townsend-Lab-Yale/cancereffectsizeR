#' Tumor allele column addition
#'
#' Several MAF files recently encountered had different representations of
#' the "tumor allele" . This function attempts to normalize all of the
#' different representations and put the tumor allele into a column called
#' "Tumor_allele".
#'
#' @param MAF MAF file with columns labels "Reference_Allele",
#'  "Tumor_Seq_Allele1", and "Tumor_Seq_Allele2"
#' @export

#TODO: add in column headers

# Different MAF files we obtained had different ways of representing
# where the "tumor allele" was located. This normalizes to a new column.

tumor_allele_adder <- function(MAF){


  #Make sure all the characters in these columns are uppercase for consistency.
  MAF$Tumor_Seq_Allele2 <- toupper(MAF$Tumor_Seq_Allele2)
  MAF$Tumor_Seq_Allele1 <- toupper(MAF$Tumor_Seq_Allele1)
  MAF$Reference_Allele <- toupper(MAF$Reference_Allele)

  #Delete rows with no information in either Tumor Seq columns.
  if(any(MAF$Tumor_Seq_Allele1=="" & MAF$Tumor_Seq_Allele2=="", na.rm=TRUE)){
    MAF <- MAF[-which(MAF$Tumor_Seq_Allele1=="" & MAF$Tumor_Seq_Allele2==""),]
  }

  #Delete rows with no information in either Tumor Seq columns.
  if(any(MAF$Tumor_Seq_Allele1=="NA" & MAF$Tumor_Seq_Allele2=="NA", na.rm=TRUE)){
    MAF <- MAF[-which(MAF$Tumor_Seq_Allele1=="NA" & MAF$Tumor_Seq_Allele2=="NA"),]
  }

  #Need to determine which column contains the tumor allele.
  tumor_allele_function <- function(MAF_input){
    ref <- MAF_input[,"Reference_Allele"]
    alt1 <- MAF_input[,"Tumor_Seq_Allele1"]
    alt2 <- MAF_input[,"Tumor_Seq_Allele2"]
    if(ref!=alt1 & ref!=alt2 & alt1==alt2){
      return(alt2) # If there wasn't an alt2, we added based on alt1, so now they are the same but not the reference.
    }
    if(ref==alt1 & ref==alt2){
      #everything is equal, return and then delete later
      return(ref)
    }
    if(ref==alt1 & ref!=alt2){
      return(alt2)
    }else{
      if(ref==alt2 & ref!=alt1){
        return(alt1)
      }else{
        if(length(which(c(alt1,alt2)==""))==1){
          if(alt1 == "" & alt2 != ref){
            return(alt2)
          }else{
            if(alt2 == "" & alt1 != ref){
              return(alt1)
            }
          }
        }
      }
    }
  }

  MAF_allele_list <- split(x = MAF[,c("Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2")],f = rep(1:nrow(MAF),each=1))
  MAF$Tumor_allele <- unlist(lapply(X = MAF_allele_list,FUN = tumor_allele_function))

  if(any(MAF$Tumor_allele == MAF$Reference_Allele, na.rm=TRUE)){
   message(paste(sum(MAF$Tumor_allele == MAF$Reference_Allele),
                 " variants have the same reference and tumor allele, \n
                 deleting from the analysis"))
    MAF[-which(MAF$Tumor_allele == MAF$Reference_Allele),]

  }
  return(MAF)
}
