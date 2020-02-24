#' Trinucleotide profiler
#'
#' Creates a dataframe with the average trinucleotide product among tumors
#' for 96 trinucleotide contexts using
#' \code{\link[deconstructSigs]{whichSignatures}}. Also outputs the individual
#' signature weights vs tumors.
#'
#' @param input.MAF MAF with mutation data
#' @param signature.choice Either "signatures.cosmic" or "signatures.nature2013".
#' @param minimum.mutation.per.tumor Set to 50, as per Rosenthal et al. (see Details)
#' @param save.figs Boolean, if T saves figures (requires \code{ggplot2})
#' to "Figures/*tumor.name*.pdf"
#' @param tumor.name If \code{save.figs = T} then you should have a character string
#' here for the save file and plot title.
#' @param sample.ID.column name of the column in \code{input.MAF} that contains unique
#' sample identifiers.
#' @param chr.column name of the column in \code{input.MAF} that contains chromosome
#' information.
#' @param pos.column name of the column in \code{input.MAF} that contains SNV position
#' information.
#' @param ref.column name of the column in \code{input.MAF} that contains reference
#' allele information.
#' @param alt.column name of the column in \code{input.MAF} that contains alternative
#' allele information.
#' @import deconstructSigs
#'
#' @details Uses the \code{deconstructSigs} package from Rosenthal, R., McGranahan, N., Herrero, J., Taylor, B. S., & Swanton, C. (2016). deconstructSigs: delineating mutational processes in single tumors distinguishes DNA repair deficiencies and patterns of carcinoma evolution. Genome Biology, 17(1), 31. https://doi.org/10.1186/s13059-016-0893-4
#' @export


# calculates the average trinucleotide context among tumors

trinuc_profile_function_with_weights <-
  function(input_MAF,
           signature_choice="signatures.cosmic",
           minimum_mutations_per_tumor=50,
           save_figs=F,tumor.name="[Tumor name here]",
           sample_ID_column="Unique_patient_identifier",
           chr_column = "Chromosome",
           pos_column = "Start_Position",
           ref_column = "Reference_Allele",
           alt_column = "Tumor_allele"){

    # ip <- as.data.frame(installed.packages()[,c(1,3:4)])
    # if (!require("deconstructSigs")) {
    #   install.packages("deconstructSigs", dependencies = TRUE)
    #   library(deconstructSigs)
    # }
    # library(deconstructSigs)

    if(length(which(colnames(input_MAF)==sample_ID_column))==0){
      warning("You require a column named \'Unique_patient_identifier\' or you need to specify the column with unique IDs for all tumors.\nReturning NULL")
      return(NULL)
    }

    input_MAF[,chr_column] <- paste("chr",trimws(input_MAF[,chr_column]),sep="")
    # should all be SNP already
    # input.MAF <- input.MAF[which(input.MAF$Variant_Type=="SNP"),]

    #We do not want the signals for selection here, so we get rid of all recurrently mutated SNP (same start position and same chromosome)
    message("Removing all recurrent mutations... ")
    # duplicated.vec.first <- duplicated(input.MAF$Start_Position)
    # duplicated.vec.last <- duplicated(input.MAF$Start_Position,fromLast=T)

    duplicated_vec_first <- duplicated(input_MAF[,c(pos_column,chr_column,alt_column)])
    duplicated_vec_last <- duplicated(input_MAF[,c(pos_column,chr_column,alt_column)],fromLast=T)

    duplicated_vec_pos <- which(duplicated_vec_first | duplicated_vec_last)
    # duplicated.vec.start <- input.MAF$Start_Position[duplicated.vec.pos]


    # MAF.test <- MAF


    # duplicated.positions <- MAF_for_analysis$Start_Position[duplicated(MAF_for_analysis$Start_Position)]

    # remove.from.eee <- NULL
    # for(eee in 1:length(duplicated.vec.pos)){
    #
    #   positions <- which(input.MAF$Start_Position==duplicated.vec.start[eee])
    #   positions <- positions[-which(positions == duplicated.vec.pos[eee])]
    #   if(!(input.MAF$Chromosome[duplicated.vec.pos[eee]] %in% input.MAF$Chromosome[positions])){
    #     remove.from.eee <- c(remove.from.eee,eee)
    #   }
    #
    # }
    # if(length(remove.from.eee)>0){
    #   duplicated.vec.pos <- duplicated.vec.pos[-remove.from.eee]
    # }
    #
    input_MAF <- input_MAF[-duplicated_vec_pos,]
    # no.recur.MAF <- MAF_for_analysis[include.vec,]

    #Finding unique tumors

    message("Finding the number of mutations per tumor")

    input_MAF <- input_MAF[which(input_MAF[, ref_column] %in% c('A', 'T', 'C', 'G') & input_MAF[, alt_column] %in% c('A', 'T', 'C', 'G')),]

    unique_patients <- unique(input_MAF[,sample_ID_column])
    tumor_mutation_number <- NULL
    for(i in 1:length(unique_patients)){
      tumor_mutation_number[i] <- nrow(input_MAF[which(input_MAF[,sample_ID_column]==unique_patients[i]),])

    }


    patients_over_mut_minimum <- unique_patients[which(tumor_mutation_number>minimum_mutations_per_tumor)]

    message(paste("Number of tumors over specified minimum mutation number of ",minimum_mutations_per_tumor,": ",length(patients_over_mut_minimum),sep=""))

    message("Cleaning input to only contain tumors above the minimum...")

    #Only keeping patients above the minimum. There has to be a much faster way for this, but this works for now...
    # keep <- NULL
    # for(i in 1:nrow(input.MAF)){
    #   if(input.MAF[i,sample.ID.column] %in% patients.over.mut.minimum){keep <- c(keep,i)}
    # }
    keep <- which(input_MAF[,sample_ID_column] %in% patients_over_mut_minimum)
    input_MAF <- input_MAF[keep,]


    message("Calculating trinucleotide mutation counts...")
    sigs_input <- deconstructSigs::mut.to.sigs.input(mut.ref = input_MAF,
                                    sample.id = sample_ID_column,
                                    chr = chr_column,
                                    pos = pos_column,
                                    ref = ref_column,
                                    alt = alt_column)

    # data("tri.counts.genome", package = "deconstructSigs")
    # data("tri.counts.exome",package = "deconstructSigs")
    # data("tri.counts.genome",package = "deconstructSigs")
    # data("signatures.cosmic",package = "deconstructSigs")
    # data("signatures.nature2013", package = "deconstructSigs")


    message("Calculating individual tumor mutational signatures...")

    if(signature_choice=="signatures.cosmic"){
      signatures_output <- list()
      for(i in 1:nrow(sigs_input)){
        signatures_output[[i]] <- deconstructSigs::whichSignatures(tumor.ref = sigs_input,
                                                  signatures.ref = signatures.cosmic,
                                                  sample.id = rownames(sigs_input)[i],
                                                  contexts.needed = TRUE,
                                                  tri.counts.method = 'exome2genome')
      }
    }
    if(signature_choice=="signatures.nature2013"){
      signatures_output <- list()
      for(i in 1:nrow(sigs_input)){
        signatures_output[[i]] <- deconstructSigs::whichSignatures(tumor.ref = sigs_input,
                                                  signatures.ref = signatures.nature2013,
                                                  sample.id = rownames(sigs_input)[i],
                                                  contexts.needed = TRUE,
                                                  tri.counts.method = 'exome2genome')
      }

    }


    # melt(rbind(signatures.output[[1]]$weights,signatures.output[[2]]$weights))

    if(length(signatures_output)>1){
      dir.create(path = "Figures", showWarnings = FALSE)
      weights_df <- signatures_output[[1]]$weights
      for(i in 2:length(signatures_output)){
        weights_df <- rbind(weights_df,signatures_output[[i]]$weights)
      }

      weights_df_melt <- reshape2::melt(weights_df)

      if(save_figs){
        viol <- ggplot2::ggplot(data=weights_df_melt,aes(x=variable,y=value)) + geom_violin() + theme_bw() +
          theme(axis.text.x = element_text(angle=90,hjust = 1)) +
          geom_jitter(width=0.25,aes(x=variable,y=value),alpha=0.2) +
          labs(x="Cosmic signature") + labs(y="Signature weight") + ggtitle(paste("Signatures in ",tumor.name,sep="")) +
          theme(plot.title = element_text(hjust = 0.5))

        # viol
        ggsave(plot = viol,filename = paste("Figures/Signature_Violin_plot_",tumor.name,".pdf",sep=""),units = "in",width = 11,height = 8)
      }



    }else{
      message("Only one signature calculated!")
      dir.create(path = "Figures", showWarnings = FALSE)
      weights_df <- signatures_output[[1]]$weights
      weights_df_melt <- reshape2::melt(weights_df)
    }




    unknown_list <- rep(NA,length(signatures_output))
    for(i in 1:length(signatures_output)){
      unknown_list[i] <- signatures_output[[i]]$unknown
    }
    message("Statistical summary of the proportion of the mutational signature in each tumor sample that is \'unknown\'")
    print(summary(unknown_list))

    averaged_product <- signatures_output[[1]]$product
    for(i in 2:length(signatures_output)){
      averaged_product <- averaged_product + signatures_output[[i]]$product
    }
    averaged_product <- averaged_product/length(signatures_output)
    averaged_product <- averaged_product/sum(averaged_product)
    # sum(averaged.product)

    ###
    #Create trinucleotide data frame
    ###
    trinuc_mutation_data <- expand.grid(mutation=c("CtoA","CtoG","CtoT","TtoA","TtoC","TtoG"),Upstream=c("A","C","G","T"),Downstream=c("A","C","G","T"))
    trinuc_mutation_data$mutated_from <- rep(NA,nrow(trinuc_mutation_data))
    trinuc_mutation_data$mutated_to <- rep(NA,nrow(trinuc_mutation_data))
    for(i in 1:nrow(trinuc_mutation_data)){
      trinuc_mutation_data$mutated_from[i] <- strsplit(as.character(trinuc_mutation_data$mutation[i]),split = "")[[1]][1]
    }
    for(i in 1:nrow(trinuc_mutation_data)){
      trinuc_mutation_data$mutated_to[i] <- strsplit(as.character(trinuc_mutation_data$mutation[i]),split = "")[[1]][4]
    }
    trinuc_mutation_data$total_count <- rep(0,nrow(trinuc_mutation_data)) ##
    trinuc_mutation_data$proportion <- rep(0,nrow(trinuc_mutation_data)) ##
    trinuc_mutation_data$section_labels <- factor(trinuc_mutation_data$mutation, labels = c("C%->%A", "C%->%G", "C%->%T","T%->%A", "T%->%C", "T%->%G"))



    for(i in 1:length(averaged_product)){
      splitname <- strsplit(colnames(averaged_product)[i],split="")[[1]]

      trinuc_mutation_data$proportion[which(trinuc_mutation_data$Upstream==splitname[1] &
                                              trinuc_mutation_data$Downstream==splitname[7] &
                                              trinuc_mutation_data$mutated_from==splitname[3] &
                                              trinuc_mutation_data$mutated_to==splitname[5])] <- averaged_product[i]
    }
    if(save_figs){
      p <- ggplot2::ggplot(data=trinuc_mutation_data, aes(Downstream, Upstream)) +
        geom_tile(aes(fill = proportion*100), colour = "white") + scale_fill_gradient(low = "white", high = "steelblue", name="Percent")
      p <- p + facet_grid(.~section_labels, labeller = label_parsed)
      p <- p +  geom_text(aes(label = round(proportion, 4)*100),size=3)
      # p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      # panel.background = element_blank(), axis.line = element_line(colour = "black"))
      p <- p + theme_bw() + theme(panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(),
                                  axis.ticks = element_blank(),
                                  strip.text=element_text(size=15),
                                  axis.title.x = element_text(size=15),
                                  axis.title.y = element_text(size=15),
                                  axis.text.x = element_text(size=12),
                                  axis.text.y=element_text(size=12),plot.title = element_text(hjust = 0.5)) +
        ggtitle(paste("Trinucleotide profile for ",tumor.name,sep=""))
      # p
      ggsave(paste("Figures/",tumor.name,"_trinuc_heatmap.pdf",sep=""),height = 2.5,width = 10)
    }

    trinuc_output <- list(trinuc_mutation_data=trinuc_mutation_data,
                          signature_weights=signatures_output)
    return(trinuc_output)
  }

