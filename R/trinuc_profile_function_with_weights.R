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
  function(input.MAF,
           signature.choice="signatures.cosmic",
           minimum.mutations.per.tumor=50,
           save.figs=F,tumor.name="[Tumor name here]",
           sample.ID.column="Unique_patient_identifier",
           chr.column = "Chromosome",
           pos.column = "Start_Position",
           ref.column = "Reference_Allele",
           alt.column = "Tumor_allele"){

    # ip <- as.data.frame(installed.packages()[,c(1,3:4)])
    # if (!require("deconstructSigs")) {
    #   install.packages("deconstructSigs", dependencies = TRUE)
    #   library(deconstructSigs)
    # }
    # library(deconstructSigs)

    if(length(which(colnames(input.MAF)==sample.ID.column))==0){
      warning("You require a column named \'Unique_patient_identifier\' or you need to specify the column with unique IDs for all tumors.\nReturning NULL")
      return(NULL)
    }

    input.MAF[,chr.column] <- paste("chr",trimws(input.MAF[,chr.column]),sep="")
    # should all be SNP already
    # input.MAF <- input.MAF[which(input.MAF$Variant_Type=="SNP"),]

    #We do not want the signals for selection here, so we get rid of all recurrently mutated SNP (same start position and same chromosome)
    message("Removing all recurrent mutations... ")
    # duplicated.vec.first <- duplicated(input.MAF$Start_Position)
    # duplicated.vec.last <- duplicated(input.MAF$Start_Position,fromLast=T)

    duplicated.vec.first <- duplicated(input.MAF[,c(pos.column,chr.column,alt.column)])
    duplicated.vec.last <- duplicated(input.MAF[,c(pos.column,chr.column,alt.column)],fromLast=T)

    duplicated.vec.pos <- which(duplicated.vec.first | duplicated.vec.last)
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
    input.MAF <- input.MAF[-duplicated.vec.pos,]
    # no.recur.MAF <- MAF_for_analysis[include.vec,]

    #Finding unique tumors

    message("Finding the number of mutations per tumor")
    unique.patients <- unique(input.MAF[,sample.ID.column])
    tumor.mutation.number <- NULL
    for(i in 1:length(unique.patients)){
      tumor.mutation.number[i] <- nrow(input.MAF[which(input.MAF[,sample.ID.column]==unique.patients[i]),])

    }


    patients.over.mut.minimum <- unique.patients[which(tumor.mutation.number>minimum.mutations.per.tumor)]

    message(paste("Number of tumors over specified minimum mutation number of ",minimum.mutations.per.tumor,": ",length(patients.over.mut.minimum),sep=""))

    message("Cleaning input to only contain tumors above the minimum...")

    #Only keeping patients above the minimum. There has to be a much faster way for this, but this works for now...
    # keep <- NULL
    # for(i in 1:nrow(input.MAF)){
    #   if(input.MAF[i,sample.ID.column] %in% patients.over.mut.minimum){keep <- c(keep,i)}
    # }
    keep <- which(input.MAF[,sample.ID.column] %in% patients.over.mut.minimum)
    input.MAF <- input.MAF[keep,]


    message("Calculating trinucleotide mutation counts...")
    sigs.input <- deconstructSigs::mut.to.sigs.input(mut.ref = input.MAF,
                                    sample.id = sample.ID.column,
                                    chr = chr.column,
                                    pos = pos.column,
                                    ref = ref.column,
                                    alt = alt.column)

    # data("tri.counts.genome", package = "deconstructSigs")
    # data("tri.counts.exome",package = "deconstructSigs")
    # data("tri.counts.genome",package = "deconstructSigs")
    # data("signatures.cosmic",package = "deconstructSigs")
    # data("signatures.nature2013", package = "deconstructSigs")


    message("Calculating individual tumor mutational signatures...")

    if(signature.choice=="signatures.cosmic"){
      signatures.output <- list()
      for(i in 1:nrow(sigs.input)){
        signatures.output[[i]] <- deconstructSigs::whichSignatures(tumor.ref = sigs.input,
                                                  signatures.ref = signatures.cosmic,
                                                  sample.id = rownames(sigs.input)[i],
                                                  contexts.needed = TRUE,
                                                  tri.counts.method = 'exome2genome')
      }
    }
    if(signature.choice=="signatures.nature2013"){
      signatures.output <- list()
      for(i in 1:nrow(sigs.input)){
        signatures.output[[i]] <- deconstructSigs::whichSignatures(tumor.ref = sigs.input,
                                                  signatures.ref = signatures.nature2013,
                                                  sample.id = rownames(sigs.input)[i],
                                                  contexts.needed = TRUE,
                                                  tri.counts.method = 'exome2genome')
      }

    }


    # melt(rbind(signatures.output[[1]]$weights,signatures.output[[2]]$weights))

    if(length(signatures.output)>1){
      dir.create(path = "Figures", showWarnings = FALSE)
      weights.df <- signatures.output[[1]]$weights
      for(i in 2:length(signatures.output)){
        weights.df <- rbind(weights.df,signatures.output[[i]]$weights)
      }

      weights.df.melt <- reshape2::melt(weights.df)

      if(save.figs){
        viol <- ggplot2::ggplot(data=weights.df.melt,aes(x=variable,y=value)) + geom_violin() + theme_bw() +
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
      weights.df <- signatures.output[[1]]$weights
      weights.df.melt <- reshape2::melt(weights.df)
    }




    unknown.list <- rep(NA,length(signatures.output))
    for(i in 1:length(signatures.output)){
      unknown.list[i] <- signatures.output[[i]]$unknown
    }
    message("Statistical summary of the proportion of the mutational signature in each tumor sample that is \'unknown\'")
    print(summary(unknown.list))

    averaged.product <- signatures.output[[1]]$product
    for(i in 2:length(signatures.output)){
      averaged.product <- averaged.product + signatures.output[[i]]$product
    }
    averaged.product <- averaged.product/length(signatures.output)
    averaged.product <- averaged.product/sum(averaged.product)
    # sum(averaged.product)

    ###
    #Create trinucleotide data frame
    ###
    trinuc.mutation_data <- expand.grid(mutation=c("CtoA","CtoG","CtoT","TtoA","TtoC","TtoG"),Upstream=c("A","C","G","T"),Downstream=c("A","C","G","T"))
    trinuc.mutation_data$mutated_from <- rep(NA,nrow(trinuc.mutation_data))
    trinuc.mutation_data$mutated_to <- rep(NA,nrow(trinuc.mutation_data))
    for(i in 1:nrow(trinuc.mutation_data)){
      trinuc.mutation_data$mutated_from[i] <- strsplit(as.character(trinuc.mutation_data$mutation[i]),split = "")[[1]][1]
    }
    for(i in 1:nrow(trinuc.mutation_data)){
      trinuc.mutation_data$mutated_to[i] <- strsplit(as.character(trinuc.mutation_data$mutation[i]),split = "")[[1]][4]
    }
    trinuc.mutation_data$total_count <- rep(0,nrow(trinuc.mutation_data)) ##
    trinuc.mutation_data$proportion <- rep(0,nrow(trinuc.mutation_data)) ##
    trinuc.mutation_data$section_labels <- factor(trinuc.mutation_data$mutation, labels = c("C%->%A", "C%->%G", "C%->%T","T%->%A", "T%->%C", "T%->%G"))



    for(i in 1:length(averaged.product)){
      splitname <- strsplit(colnames(averaged.product)[i],split="")[[1]]

      trinuc.mutation_data$proportion[which(trinuc.mutation_data$Upstream==splitname[1] &
                                              trinuc.mutation_data$Downstream==splitname[7] &
                                              trinuc.mutation_data$mutated_from==splitname[3] &
                                              trinuc.mutation_data$mutated_to==splitname[5])] <- averaged.product[i]
    }
    if(save.figs){
      p <- ggplot2::ggplot(data=trinuc.mutation_data, aes(Downstream, Upstream)) +
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

    trinuc.output <- list(trinuc.mutation_data=trinuc.mutation_data,
                          signature.weights=signatures.output)
    return(trinuc.output)
  }

