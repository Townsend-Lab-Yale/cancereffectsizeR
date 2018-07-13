#' Unique tumor addition function
#'
#' Adds a column called "Unique_patient_identifier" to your MAF file,
#' useful if TCGA codes have more data than sample and you want summary
#' statistics of the number of mutations per tumor.
#'
#' @param MAF.file MAF file with "Tumor_Sample_Barcode" column
#' @param non.TCGA.characters.to.keep numeric number of characters to keep
#' in "Tumor_Sample_Barcode" column for sample names that
#' do not start with "TCGA". If non-numeric object then defaults to
#' all characters
#' @param sum.stats Boolean, T prints summary() of number of mutations per tumor
#' @param figures Boolean, T saves histogram of sum.stats and saves as "Figures/mutations_per_tumor.pdf"
#' @export



# Adds column to the MAF that contains the unique tumor ID from the $Tumor_Sample_Barcode

unique_tumor_addition_function <- function(MAF.file,non.TCGA.characters.to.keep="all",sum.stats=T,figures=F){
  this.maf <- MAF.file
  TCGA.tumors <- grep(pattern = "TCGA",x = this.maf$Tumor_Sample_Barcode)
  other.tumors <- setdiff(1:nrow(this.maf),TCGA.tumors)

  first.12 <- function(string_to_12){
    return(paste(unlist(strsplit(string_to_12,split = ""))[1:12],collapse = ""))
  }

  if(class(non.TCGA.characters.to.keep)=="numeric"){
    first.other <- function(string_to_other){
      return(paste(unlist(strsplit(string_to_other,split = ""))[1:non.TCGA.characters.to.keep],collapse = ""))
    }

    this.maf$Unique_patient_identifier <- NA

    this.maf$Unique_patient_identifier[TCGA.tumors] <- unlist(lapply(this.maf$Tumor_Sample_Barcode[TCGA.tumors],first.12))
    this.maf$Unique_patient_identifier[other.tumors] <- unlist(lapply(this.maf$Tumor_Sample_Barcode[other.tumors],first.other))


  }else{
    this.maf$Unique_patient_identifier <- NA

    this.maf$Unique_patient_identifier[TCGA.tumors] <- unlist(lapply(this.maf$Tumor_Sample_Barcode[TCGA.tumors],first.12))
    this.maf$Unique_patient_identifier[other.tumors] <- this.maf$Tumor_Sample_Barcode[other.tumors]

  }

  if(sum.stats){
    unique.patients <- unique(this.maf$Unique_patient_identifier)

    tumor.mutation.number <- NULL
    for(i in 1:length(unique.patients)){
      tumor.mutation.number[i] <- nrow(this.maf[which(this.maf$Unique_patient_identifier==unique.patients[i]),])

    }

    if(figures){
      pdf(file = "Figures/mutations_per_tumor.pdf")
      hist(tumor.mutation.number,breaks=100,xlab="Number of mutations per tumor",main="Histogram of the number of mutations in each tumor")
      dev.off()
    }

    message("Summary statistics of the number of mutations per unique tumor:")
    print(summary(tumor.mutation.number))

  }

  return(this.maf)
}
