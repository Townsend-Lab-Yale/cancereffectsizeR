#' Selection intensity calculation (SNV)
#'
#' Function that calculates the selection intensity of single nucleotide
#' variants in cancer exome data.
#'
#' @param genes_for_analysis String of genes, e.g. c("KRAS","EGFR","CTNNB1"), to calculate effect size.
#' If "all", the function calculates the effect size of all genes within the
#' data frame.
#' @param MAF_for_analysis Data frame in MAF format
#' (\url{https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/})
#' with single nucleotide substitution data. Must contain columns with
#' gene names, chromosome, chromosomal position, reference allele, and mutant allele.
#' @param this.substitution Optional parameter to introduce a custom substitution
#' into the mutation rate calculation, e.g. like KRAS G12C here \url{https://www.nature.com/articles/s41388-017-0105-z}.
#' @param trinuc.mutation_data Data frame with trinucleotide mutation proportions for the tumor type.
#' @param mut_rates Vector with gene-level mutation rate data.
#' @param dndscv_siggenes dndscv data to be integrated into results data table.
#' @param tumor.number Number of tumor samples.
#' Useful if you split the analysis up into several independent runs. If \code{NULL}
#' then automatically sets to number of unique samples in \code{MAF_for_analysis}.
#' @param output_from_mainMAF the output from a previous selection run with all the
#' mutation and selection data, useful for bootstrapping because it speeds up the
#' analysis (no longer need to look up mutational context, etc.)
#' @param sample_ID_column column in MAF with sample ID data
#' @param gene_ID_column column in MAF with gene ID data
#' @param ref_column column in MAF with reference allele data
#' @param alt_column column in MAF with alternative allele data
#' @param pos_column column in MAF with chromosome nucleotide location data
#' @param chr_column column in MAF with chromosome data
#'
#'
#' @import dndscv
#' @import BSgenome
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @export


selection_intensity_calculation <- function(genes_for_analysis="all",
                                            MAF_for_analysis,
                                            this.substitution=c("fake_gene",41,"T"),
                                            trinuc.mutation_data,
                                            mut_rates,
                                            dndscv_siggenes,
                                            tumor.number=NULL,
                                            trinuc_count_matrix=NULL,
                                            output_from_mainMAF=NULL,
                                            sample_ID_column="Unique_patient_identifier",
                                            gene_ID_column="Gene_name",
                                            ref_column="Reference_Allele",
                                            alt_column="Tumor_allele",
                                            chr_column="Chromosome",
                                            pos_column="Start_Position",
                                            every_gene_verbose=T,
                                            ref_cds_object=NULL,custom_BSgenome_selection=NULL){

  # data("all_gene_trinuc_data",package = "cancereffectsizeR")
  data("AA_translations", package = "cancereffectsizeR")



  if(genes_for_analysis=="all"){
    genes_for_analysis <- unique(MAF_for_analysis[,gene_ID_column])
  }

  MAF_for_analysis = MAF_for_analysis[which(
    (MAF_for_analysis[,ref_column] %in% c("A","T","G","C")) &
      (MAF_for_analysis[,alt_column] %in% c("A","T","G","C"))),]

  flip.function <- function(nucleotide){
    if(nucleotide=="A"){
      return("T")
    }
    if(nucleotide=="T"){
      return("A")
    }
    if(nucleotide=="C"){
      return("G")
    }
    if(nucleotide=="G"){
      return("C")
    }
    if(nucleotide=="N"){
      return("N")
    }
  }


  lambda.calc <- function(n.one,n.zero){
    return(-log(n.zero/(n.one+n.zero))) # maximum likelihood of 1 or more substitutions at a site
  }



  if(is.null(tumor.number)){
    warning("You did not provide the number of tumors in this dataset!\nUsing the unique tumors in the \'Unique_patient_identifier\' column of the MAF_for_analysis input.\nIf you broke the MAF into separate files, you will be missing the true number of tumors in this dataset!")
    tumor.number <- length(unique(MAF_for_analysis[,sample_ID_column]))
  }else{
    tumor.number <- tumor.number
  }

  # data("refcds_hg19", package="dndscv")
  RefCDS <- ref_cds_object
  # RefCDS[[1]]
  # gr_genes
  gene_list  <- sapply(RefCDS, function(x) x$gene_name)


  #
  # Loading data into the environment now so it is only called once per function.

  if(length(which(MAF_for_analysis[,gene_ID_column]=="WI2-3308P17.2"))>0){
    MAF_for_analysis <- MAF_for_analysis[-which(MAF_for_analysis[,gene_ID_column]=="WI2-3308P17.2"),]
  }



  MAF_for_analysis$End_Position <- MAF_for_analysis[,pos_column]

  if(length(which(colnames(MAF_for_analysis)=="Start_position"))>0){
    colnames(MAF_for_analysis)[which(colnames(MAF_for_analysis)=="Start_position")] <- "Start_Position"
  }
  AA_translations <- as.matrix(AA_translations)
  rownames(AA_translations) <- AA_translations[,"Nucs"]

  # Creating two class-specific matrices of the trinucleotide context because it saves a lot of time looping through later
  trinuc.mutation_data.char <- as.matrix(trinuc.mutation_data[,c("mutation","Upstream","Downstream","mutated_from","mutated_to","section_labels")])
  trinuc.mutation_data.num <- as.matrix(trinuc.mutation_data[,c("total_count","proportion")])

  namesvec <- NULL
  for(i in 1:nrow(trinuc.mutation_data)){
    namesvec[i] <- paste0(trinuc.mutation_data$mutated_from[i],trinuc.mutation_data$Upstream[i],trinuc.mutation_data$mutated_to[i],trinuc.mutation_data$Downstream[i],collapse = "")
  }
  # namesvec
  rownames(trinuc.mutation_data.num) <- namesvec




  all.muts <- as.data.frame(matrix(data=NA,nrow=0,ncol=16),header=T) #Data frame to hold mutation data
  colnames(all.muts) <- c("Gene",
                          "AA_Pos",
                          "Nucleotide_position",
                          "Nuc_Ref",
                          "Nuc_Change",
                          "AA_Ref",
                          "AA_Change",
                          "gamma",
                          "selection_intensity",
                          "prevalence_among_tumors",
                          "mutation_rate",
                          "gene_AA_size",
                          "dndscv_p",
                          "dndscv_q",
                          "Prop_tumors_with_specific_mut",
                          "Prop_of_tumors_with_this_gene_mutated") #Column names for this data frame


  #Create a data frame to hold all the mutational data.
  mutation.data <- as.data.frame(matrix(data=NA,nrow=0,ncol=30),header=T)
  colnames(mutation.data) <- c("Gene",
                               "Gene_isoform",
                               "Gene_size",
                               "Percent_A",
                               "Percent_T",
                               "Percent_G",
                               "Percent_C",
                               "Nucleotide_Gene_Pos",
                               "Nucleotide_chromosome_position",
                               "Chromosome",
                               "Reference_Nucleotide",
                               "Alternative_Nucleotide",
                               "Reference_Count",
                               "Alternative_Count",
                               "Tumor_origin",
                               "Nucleotide_change_tally",
                               "Nucleotide_mutation_rate",
                               "Amino_acid_position",
                               "Amino_acid_codon",
                               "Codon_position",
                               "Amino_acid_reference",
                               "Amino_acid_alternative",
                               "Amino_acid_mutation_rate",
                               "Amino_acid_change_tally",
                               "Gamma",
                               "MAF_location",
                               "Nucleotide_trinuc_context",
                               "gene_level_synonymous_mutation_rate",
                               "strand",
                               "trinucs")

  #matrix to hold trinucleotide count data
  if(is.null(trinuc_count_matrix)){
    trinuc.count.matrix <- matrix(nrow=96,ncol=length(unique(MAF_for_analysis[,gene_ID_column])),data=NA)
    colnames(trinuc.count.matrix) <- unique(MAF_for_analysis[,gene_ID_column])
    rownames(trinuc.count.matrix) <- namesvec
  }else{trinuc.count.matrix <- trinuc_count_matrix}

  #Going to cut all substitutions out of the analysis that do not match the reference genome. Will keep track.
  if(is.null(output_from_mainMAF)){
    mutations.to.cut <- as.data.frame(matrix(data=NA,nrow=1,ncol=3))
    colnames(mutations.to.cut) <- c(gene_ID_column,chr_column,pos_column)
  }

  ####Loop to iterate through every gene----
  for(zzz in 1:length(genes_for_analysis)){


    this.gene <- genes_for_analysis[zzz] #The specific gene for this run is... <-
    if(every_gene_verbose){print(paste("This gene:",this.gene))}
    #If the name of this gene is "Unknown", tell us!
    if(this.gene=="Unknown" | this.gene=="Indeterminate"){
      print("This gene is named UNKNOWN or Indeterminate!")
      next
    }



    # assume MAF is just SNV at this point but you can uncomment this block
    # if not.
    this.gene_MAF <- MAF_for_analysis[which(MAF_for_analysis[,gene_ID_column]==this.gene #& MAF_for_analysis$Variant_Type=="SNP"
                                            #                                         # (
                                            #                                         #   MAF_for_analysis$Variant_Classification=="Nonsense_Mutation" |
                                            #                                         #     MAF_for_analysis$Variant_Classification=="Missense_Mutation" |
                                            #                                         #     MAF_for_analysis$Variant_Classification=="Silent" |
                                            #                                         #     MAF_for_analysis$Variant_Classification=="Splice_Site" |
                                            #                                         #     MAF_for_analysis$Variant_Classification=="Nonstop_Mutation")
    ),] #get a gene-specific MAF file with just the SNPs we want to analyze




    #Make sure this gene is in the MAF file, and then...
    if(nrow(this.gene_MAF)>0){


      if(is.null(output_from_mainMAF)){


        this.strand <- RefCDS[[which(gene_list==this.gene)]]$strand
        if(this.strand==-1){this.strand <- "-"}
        if(this.strand==1){this.strand <- "+"}


        seq.exon.start.vec <- RefCDS[[which(gene_list==this.gene)]]$intervals_cds[,1] - 1
        seq.exon.end.vec <- RefCDS[[which(gene_list==this.gene)]]$intervals_cds[,2]

        ###Getting isoform sequence data ----
        if(this.strand=="+"){
          myseq <- paste(as.character(BSgenome::getSeq(custom_BSgenome_selection,paste("chr",RefCDS[[which(gene_list==this.gene)]]$chr,sep=""),start=seq.exon.start.vec+1,end=seq.exon.end.vec,strand="+")),collapse = "")


          #If this isoform strand is "-", we need to flip the positions and strand to make it "+"
        }else{

          myseq <- paste(rev(as.character(BSgenome::getSeq(custom_BSgenome_selection,paste("chr",RefCDS[[which(gene_list==this.gene)]]$chr,sep=""),start=seq.exon.start.vec+1,end=seq.exon.end.vec,strand="-"))),collapse = "")


        }
        ###replace it with the mutation----

        if(this.substitution[1]==this.gene){
          tmp.string <- strsplit(myseq,"")[[1]]
          tmp.string[as.numeric(this.substitution[2])] <- this.substitution[3]
          myseq <- paste(tmp.string,collapse = "")
        }





        tmp.string <- strsplit(myseq,"")[[1]] #temporary string with our sequence split up


        myseq.codons <- paste0(tmp.string[c(T,F,F)],tmp.string[c(F,T,F)],tmp.string[c(F,F,T)]) #every 3 as a group for the codon designations


        AA.code <- rep(NA,length=length(myseq.codons)) #The amino acid code will go here

        #Loop to give the amino acid designations for each codon. Can make this a function to speed it up. TODO
        for(i in 1:length(myseq.codons)){
          AA.code[i] <- AA_translations[which(AA_translations[,1]==myseq.codons[i]),4]
        }

        #AA.code.collapsed <- paste(AA.code,collapse = "") #collapse the amino acid code #not used anymore

        #If there is no mutation rate data for this gene, tell us.
        if(length(which(names(mut_rates)==this.gene))==0){
          print(paste("There is no mutation rate data for gene ",this.gene))
          next
        }

        #Get the synonymous mutation rate from the dndscv output.

        # gene_level_synonymous_mutation_rate <- mut_rates$r_x_X[which(mut_rates$gene==this.gene)]
        gene_level_synonymous_mutation_rate <- mut_rates[this.gene]


        #What proportion of all the mutation rates is each specific mutation
        # mut.proportion <- mut.rates/sum(mut.rates)


        ##turn the sequence into a vector
        myseq.split <- unlist(strsplit(myseq,split=""))


        #Chromosome for this gene
        # chromosome <- LabReference$chrom[this.gene.ref.loc.choice]
        chromosome <-paste("chr",RefCDS[[which(gene_list==this.gene)]]$chr,sep="")

        #Sometimes the same gene name is listed for multiple chromosomes. This throws off the sequencing, reduce the MAF for this gene to the gene chosen by isoform search
        this.gene_MAF <- this.gene_MAF[which(this.gene_MAF[,chr_column]==paste(unlist(strsplit(chromosome,split=""))[4:length(unlist(strsplit(chromosome,split="")))],collapse ="")),]


        # If there are no more mutations in the MAF, need to skip to the next gene.
        # TODO: instread of relying on previous processing for Gene_name, we should figure out ourselves given Chrom + Start_pos data
        if(nrow(this.gene_MAF)==0){
          next
        }

        # cutting out substitutions that do not match the reference genome
        matcher.vec <- as.vector(this.gene_MAF[,ref_column] == as.character(BSgenome::getSeq(custom_BSgenome_selection,chromosome,start=this.gene_MAF[,pos_column],end=this.gene_MAF$End_Position,strand="+")))

        if(length(which(matcher.vec==FALSE))>0){
          not.matching <- which(matcher.vec==F)

          mutations.to.cut <- rbind(mutations.to.cut,this.gene_MAF[not.matching,c(gene_ID_column,chr_column,pos_column)])
          this.gene_MAF <- this.gene_MAF[-not.matching,]
          if(nrow(this.gene_MAF)==0){
            print(paste("After removing inconsistencies between the reference genome and the MAF file, there are no substitutions left to analyze."))
            next
          }
        }


        #vectors to hold the position data at the gene-level and chromosome-level
        gene.pos <- rep(NA,length=((seq.exon.end.vec[length(seq.exon.end.vec)]-seq.exon.start.vec[1])+2))
        chrom.gene.pos <- rep(NA,length=((seq.exon.end.vec[length(seq.exon.end.vec)]-seq.exon.start.vec[1])+2))

        #Vector of the total sequence from start to end, including introns
        if(this.strand=="-"){
          pos.sequence  <- rev(strsplit(as.character(BSgenome::getSeq(custom_BSgenome_selection,chromosome,start=seq.exon.start.vec[1],end=seq.exon.end.vec[length(seq.exon.end.vec)]+1,strand="-")),split="")[[1]])
        }else{
          pos.sequence  <- strsplit(as.character(BSgenome::getSeq(custom_BSgenome_selection,chromosome,start=seq.exon.start.vec[1],end=seq.exon.end.vec[length(seq.exon.end.vec)]+1,strand="+")),split="")[[1]]
        }

        #All the chromosome positions in the gene, including introns
        chrom.pos <- seq.exon.start.vec[1]:(seq.exon.end.vec[length(seq.exon.end.vec)]+1)

        #All the positions in the chromosome containing exons from this gene
        nuc.pos.vec <- c()
        for(i in 1:length(seq.exon.start.vec)){
          nuc.pos.vec <- c(nuc.pos.vec,(seq.exon.start.vec[i]+1):seq.exon.end.vec[i])
        }

        # position.matrix[,1] <- nuc.pos.vec #nucleotide position within the chromosome

        # positions within the chromosome vector containing the gene vector
        these.positions <- which(chrom.pos %in% nuc.pos.vec)

        #fill in the gene-level position information, based on its position in the chromosome
        for(i in 1:length(nuc.pos.vec)){  #This could be made a function to speed it up #TODO
          # this.pos <- which(chrom.pos==nuc.pos.vec[i])
          this.pos <- these.positions[i]
          chrom.gene.pos[this.pos] <- nuc.pos.vec[i]
          if(this.strand=="+"){
            gene.pos[this.pos] <- i
          }else{
            gene.pos[this.pos] <- (length(nuc.pos.vec)+1)-i
          }
        }


        ##If the substitution is here, replace to determine new trinucleotide context.
        if(this.substitution[1]==this.gene){
          pos.sequence[which(gene.pos==this.substitution[2])] <- this.substitution[3]
        }


        #Generate the trinucleotide context at every position. The NA at the end are placeholders and are not at postions used in calculations
        minus.vec <- c(pos.sequence[-1],NA)
        plus.vec <- c(NA,pos.sequence[-length(pos.sequence)])
        if(this.strand=="-"){
          trinuc.pos <- paste(minus.vec,pos.sequence,plus.vec,sep="")
        }else{
          trinuc.pos <- paste(plus.vec,pos.sequence,minus.vec,sep="")
        }

        if("N" %in% c(minus.vec,pos.sequence,plus.vec)){
         paste(this.gene, "contains 'N' within the nucleotides of the gene, according to positions specified in the RefCDS object and the customBSgenome sequence. Skipping this gene")
          next
        }

        ##Nucleotide matrix----
        #Constructing a matrix to hold all the mutation prevalence data
        mut.matrix <- matrix(nrow=4,
                             ncol=length(myseq.split),
                             data=NA)
        rownames(mut.matrix) <- c("A","T","G","C") #Rows will be the possible nucleotides
        colnames(mut.matrix) <- paste("Pos. ",seq(1,length(myseq.split),by=1)," Ref: ",myseq.split,sep="") #columns will be positions in the gene

        mut.matrix.trinuc <- matrix(nrow=nrow(mut.matrix),
                                    ncol=ncol(mut.matrix),
                                    data=NA)
        rownames(mut.matrix.trinuc) <- c("A","T","G","C") #Rows will be the possible nucleotides
        colnames(mut.matrix.trinuc) <- paste("Pos. ",seq(1,length(myseq.split),by=1)," Ref: ",myseq.split,sep="") #columns will be positions in the gene



        #Matrix that holds the mutation data for same-nucleotide mutations. This has zero observed rate
        zeros.mat <- as.data.frame(matrix(nrow=1,ncol=4,data=0))
        colnames(zeros.mat) <- c("AtoA","TtoT","GtoG","CtoC")
        # mut.proportion.combined <-  merge(mut.proportion.combined,zeros.mat) #add this to the previous mutation rate data frame

        #Fill in the matrix with all the possible mutation proportions
        position.vec <- which(!(is.na(gene.pos)))
        if(this.strand=="-"){position.vec <- position.vec[order(-position.vec)]}

        trinuc.mutation_data.num.this.gene <- trinuc.mutation_data.num
        # trinuc.data.this.gene <- trinuc.mutation_data
        # trinuc.data.this.gene$total_count <- 0
        trinuc.mutation_data.num.this.gene[,"total_count"] <- 0

        #Loop to count the trinucleotide contexts of this gene, and assign the site-specific mutation rate
        for(i in 1:length(myseq.split)){   # I feel like this could be made into a function and sped up with the lapply function. TODO
          # this.pos <- which(gene.pos==i)


          this.trinuc.split <- unlist(strsplit(trinuc.pos[position.vec[i]],split = ""))

          if(this.trinuc.split[2]=="C" | this.trinuc.split[2]=="T"){

            #counting the number of times this trinuc was mutated in this



            mut.matrix[this.trinuc.split[2],i] <- 0
            if(this.trinuc.split[2]=="C"){
              trinuc.mutation_data.num.this.gene[c(paste0(this.trinuc.split[2],this.trinuc.split[1],"A",this.trinuc.split[3],collapse = ""),
                                                   paste0(this.trinuc.split[2],this.trinuc.split[1],"G",this.trinuc.split[3],collapse = ""),
                                                   paste0(this.trinuc.split[2],this.trinuc.split[1],"T",this.trinuc.split[3],collapse = "")),
                                                 "total_count"] <-           trinuc.mutation_data.num.this.gene[c(paste0(this.trinuc.split[2],this.trinuc.split[1],"A",this.trinuc.split[3],collapse = ""),
                                                                                                                  paste0(this.trinuc.split[2],this.trinuc.split[1],"G",this.trinuc.split[3],collapse = ""),
                                                                                                                  paste0(this.trinuc.split[2],this.trinuc.split[1],"T",this.trinuc.split[3],collapse = "")),
                                                                                                                "total_count"]+1
              mut.matrix[c("A","G","T"),i] <- c(trinuc.mutation_data.num[c(paste0("C",this.trinuc.split[1],"A",this.trinuc.split[3],collapse = ""),
                                                                           paste0("C",this.trinuc.split[1],"G",this.trinuc.split[3],collapse = ""),
                                                                           paste0("C",this.trinuc.split[1],"T",this.trinuc.split[3],collapse = "")),
                                                                         "proportion"])
              mut.matrix.trinuc[c("A","G","T"),i] <- c(paste0("C",this.trinuc.split[1],"A",this.trinuc.split[3],collapse = ""),
                                                       paste0("C",this.trinuc.split[1],"G",this.trinuc.split[3],collapse = ""),
                                                       paste0("C",this.trinuc.split[1],"T",this.trinuc.split[3],collapse = ""))

            }else{
              trinuc.mutation_data.num.this.gene[c(paste0(this.trinuc.split[2],this.trinuc.split[1],"A",this.trinuc.split[3],collapse = ""),
                                                   paste0(this.trinuc.split[2],this.trinuc.split[1],"C",this.trinuc.split[3],collapse = ""),
                                                   paste0(this.trinuc.split[2],this.trinuc.split[1],"G",this.trinuc.split[3],collapse = "")),
                                                 "total_count"] <-           trinuc.mutation_data.num.this.gene[c(paste0(this.trinuc.split[2],this.trinuc.split[1],"A",this.trinuc.split[3],collapse = ""),
                                                                                                                  paste0(this.trinuc.split[2],this.trinuc.split[1],"C",this.trinuc.split[3],collapse = ""),
                                                                                                                  paste0(this.trinuc.split[2],this.trinuc.split[1],"G",this.trinuc.split[3],collapse = "")),"total_count"]+1

              # mut.matrix[c("A","C","G"),i] <- tri.nuc.flip.total.prop[positions.in.barplot]
              mut.matrix[c("A","C","G"),i] <- c(trinuc.mutation_data.num[c(paste0("T",this.trinuc.split[1],"A",this.trinuc.split[3],collapse = ""),
                                                                           paste0("T",this.trinuc.split[1],"C",this.trinuc.split[3],collapse = ""),
                                                                           paste0("T",this.trinuc.split[1],"G",this.trinuc.split[3],collapse = "")),
                                                                         "proportion"])

              mut.matrix.trinuc[c("A","C","G"),i] <- c(paste0("T",this.trinuc.split[1],"A",this.trinuc.split[3],collapse = ""),
                                                       paste0("T",this.trinuc.split[1],"C",this.trinuc.split[3],collapse = ""),
                                                       paste0("T",this.trinuc.split[1],"G",this.trinuc.split[3],collapse = ""))


            }
          }else{



            mut.matrix[this.trinuc.split[2],i] <- 0

            if(this.trinuc.split[2]=="G"){
              trinuc.mutation_data.num.this.gene[c(paste0("C",flip.function(this.trinuc.split[3]),"A",flip.function(this.trinuc.split[1]),collapse = ""),
                                                   paste0("C",flip.function(this.trinuc.split[3]),"G",flip.function(this.trinuc.split[1]),collapse = ""),
                                                   paste0("C",flip.function(this.trinuc.split[3]),"T",flip.function(this.trinuc.split[1]),collapse = "")),"total_count"] <- trinuc.mutation_data.num.this.gene[c(paste0("C",flip.function(this.trinuc.split[3]),"A",flip.function(this.trinuc.split[1]),collapse = ""),
                                                                                                                                                                                                                 paste0("C",flip.function(this.trinuc.split[3]),"G",flip.function(this.trinuc.split[1]),collapse = ""),
                                                                                                                                                                                                                 paste0("C",flip.function(this.trinuc.split[3]),"T",flip.function(this.trinuc.split[1]),collapse = "")),"total_count"]+1

              mut.matrix[c("T","C","A"),i] <- c(trinuc.mutation_data.num[c(paste0("C",flip.function(this.trinuc.split[3]),"A",flip.function(this.trinuc.split[1]),collapse = ""),
                                                                           paste0("C",flip.function(this.trinuc.split[3]),"G",flip.function(this.trinuc.split[1]),collapse = ""),
                                                                           paste0("C",flip.function(this.trinuc.split[3]),"T",flip.function(this.trinuc.split[1]),collapse = "")),"proportion"])

              mut.matrix.trinuc[c("T","C","A"),i] <- c(paste0("C",flip.function(this.trinuc.split[3]),"A",flip.function(this.trinuc.split[1]),collapse = ""),
                                                       paste0("C",flip.function(this.trinuc.split[3]),"G",flip.function(this.trinuc.split[1]),collapse = ""),
                                                       paste0("C",flip.function(this.trinuc.split[3]),"T",flip.function(this.trinuc.split[1]),collapse = ""))

            }else{

              trinuc.mutation_data.num.this.gene[c(paste0("T",flip.function(this.trinuc.split[3]),"A",flip.function(this.trinuc.split[1]),collapse = ""),
                                                   paste0("T",flip.function(this.trinuc.split[3]),"C",flip.function(this.trinuc.split[1]),collapse = ""),
                                                   paste0("T",flip.function(this.trinuc.split[3]),"G",flip.function(this.trinuc.split[1]),collapse = "")),"total_count"] <- trinuc.mutation_data.num.this.gene[c(paste0("T",flip.function(this.trinuc.split[3]),"A",flip.function(this.trinuc.split[1]),collapse = ""),
                                                                                                                                                                                                                 paste0("T",flip.function(this.trinuc.split[3]),"C",flip.function(this.trinuc.split[1]),collapse = ""),
                                                                                                                                                                                                                 paste0("T",flip.function(this.trinuc.split[3]),"G",flip.function(this.trinuc.split[1]),collapse = "")),"total_count"]+1

              mut.matrix[c("T","G","C"),i] <- c(trinuc.mutation_data.num[c(paste0("T",flip.function(this.trinuc.split[3]),"A",flip.function(this.trinuc.split[1]),collapse = ""),
                                                                           paste0("T",flip.function(this.trinuc.split[3]),"C",flip.function(this.trinuc.split[1]),collapse = ""),
                                                                           paste0("T",flip.function(this.trinuc.split[3]),"G",flip.function(this.trinuc.split[1]),collapse = "")),
                                                                         "proportion"])

              mut.matrix.trinuc[c("T","G","C"),i] <- c(paste0("T",flip.function(this.trinuc.split[3]),"A",flip.function(this.trinuc.split[1]),collapse = ""),
                                                       paste0("T",flip.function(this.trinuc.split[3]),"C",flip.function(this.trinuc.split[1]),collapse = ""),
                                                       paste0("T",flip.function(this.trinuc.split[3]),"G",flip.function(this.trinuc.split[1]),collapse = ""))
            }
          }
        }

        trinuc.count.matrix[,this.gene] <- trinuc.mutation_data.num.this.gene[,"total_count"]

        #trinuc.data.this.gene$mut_counts_times_proportion <- trinuc.data.this.gene$total_count*trinuc.data.this.gene$proportion
        # Do not need to generate the gene in the future if we calculate the mean rate in the gene once.
        # This will speed up the bootstrap analysis later.
        to.mean <-NULL
        for(w in 1:32){
          # to.mean[w] <- trinuc.data.this.gene$total_count[w*3] * sum(trinuc.data.this.gene$proportion[c(w*3,(w*3)-1,(w*3)-2)])
          to.mean[w] <- trinuc.mutation_data.num.this.gene[w*3,"total_count"] * sum(trinuc.mutation_data.num.this.gene[c(w*3,(w*3)-1,(w*3)-2),"proportion"])
        }
        mean.nuc.rate.for.this.gene <- (sum(to.mean)/(sum(trinuc.mutation_data.num.this.gene[,"total_count"])/3))

        #Now that we have the mean, we can calculate the mutation rate of nucleotides without generating the gene.

        #Find the sum of all rates at every position
        nuc.mutation.rate.vec <- rep(NA,length=length(myseq.split))
        for(i in 1:length(myseq.split)){
          nuc.mutation.rate.vec[i] <- sum(mut.matrix[,i])
        }

        #matrix that has all the synonymous mutation rates
        normalized.mut.matrix <-  (mut.matrix/mean(nuc.mutation.rate.vec))*gene_level_synonymous_mutation_rate

        ##Amino acid matrix----
        #Matrix that stores rates for each amino acid mutation
        AA.mut.matrix <- matrix(nrow=21,
                                ncol=length(AA.code),
                                data=0)
        rownames(AA.mut.matrix) <- unique(AA_translations[,3]) #Rows are all the possible nucleotides

        #Names of the columns of this matrix will be the amino acid at that location
        aa.code.short <- rep(NA,length=length(AA.code))
        for(i in 1:length(aa.code.short)){
          aa.code.short[i] <- AA_translations[which(AA_translations[,4]==AA.code[i])[1],3]
        }
        colnames(AA.mut.matrix) <- paste("Pos. ",seq(1,length(AA.code),by=1)," Ref: ",aa.code.short,sep="")

        AA.matrix.trinuc <- matrix(nrow=nrow(AA.mut.matrix),
                                   ncol=ncol(AA.mut.matrix),
                                   data=NA)
        rownames(AA.matrix.trinuc) <- rownames(AA.mut.matrix)

        #########
        #Probabilities of each amino acid mutation----
        #########

        #Number/positional designation for each codon in terms of nucleotide position
        codon.numbers <- rep(1:length(myseq.codons),each=3)
        nuc.list <- c("A","T","G","C")

        #Loop to calculate mutation rates at each codon to other codons.
        for(i in 1:length(myseq.split)){
          this.codon <- myseq.codons[codon.numbers[i]] #What is this codon?
          #see which position in the codon we are looking at
          if(i%%3==1){pos <- 1};if(i%%3==2){pos <- 2};if(i%%3==0){pos <- 3} #gives back position within the codon
          this.codon.expanded <- unlist(strsplit(this.codon,split="")) #breaks this codon into individual nucleotides
          # if(this.codon.expanded[2]=="G" | this.codon.expanded[2]=="A"){
          #
          # }
          mutations <- nuc.list[which(nuc.list!=this.codon.expanded[pos])] #possible mutations at this location within the codon

          #Loop to go through all the possible mutations at this position in the codon and add up rates
          for(j in 1:3){
            this.mutation <- mutations[j] #This specific mutation is...
            this.codon.expanded[pos] <- this.mutation #Making this codon...
            this.codon.collapsed <- paste0(this.codon.expanded,collapse="") #which we collapse...
            new.aa <- AA_translations[this.codon.collapsed,3] #resulting in this new amino acid...
            AA.mut.matrix[new.aa,codon.numbers[i]] <- AA.mut.matrix[new.aa,codon.numbers[i]] + normalized.mut.matrix[this.mutation,i] #adding the rate this happens to total rates matrix
            # AA.matrix.trinuc[new.aa,codon.numbers[i]] <- ifelse(is.na(AA.matrix.trinuc[new.aa,codon.numbers[i]]),mut.matrix.trinuc[this.mutation,i],paste(AA.matrix.trinuc[new.aa,codon.numbers[i]],mut.matrix.trinuc[this.mutation,i],sep=":"))
          }
        }

        ##Matrix of the total rates of amino acid mutations


        #matrix that stores info on the complimentary strand
        strand.switch <- as.data.frame(matrix(data=c("A","T","G","C"),ncol=1,nrow=4),stringsAsFactors = FALSE)
        rownames(strand.switch) <- c("T","A","C","G")


        ##Amino acid matrix that stores gamma information
        AA.gamma.matrix <- matrix(nrow=21,
                                  ncol=length(AA.code),
                                  data=NA)

        rownames(AA.gamma.matrix) <- rownames(AA.mut.matrix)
        colnames(AA.gamma.matrix) <- colnames(AA.mut.matrix)

        ##Amino acid that stores mutation incidence data
        AA.tally.matrix <- matrix(nrow=21,
                                  ncol=length(AA.code),
                                  data=0)

        rownames(AA.tally.matrix) <- rownames(AA.mut.matrix)
        colnames(AA.tally.matrix) <- colnames(AA.mut.matrix)


        #Matrix to store the incidence of mutations at the nucleotide level
        Nuc.tally.matrix <- matrix(nrow=4,
                                   ncol=length(myseq.split),
                                   data=0)

        rownames(Nuc.tally.matrix) <- rownames(mut.matrix)
        colnames(Nuc.tally.matrix) <- colnames(mut.matrix)


        #Matrix to tally up the mutations outside of exons
        outside.nucs <- as.data.frame(matrix(nrow=nrow(this.gene_MAF),ncol=5,data=NA))
        colnames(outside.nucs) <- c("nuc_pos","ref","change","tally","mutation_rate")
        outside.nucs$tally <- 0

        ##Setting up the data frame to store (and eventually compile) all the analysis
        mutation.data.to.add <- as.data.frame(matrix(data=NA,nrow=nrow(this.gene_MAF),ncol=31),header=T)
        colnames(mutation.data.to.add) <- c("Gene",
                                            "Gene_isoform",
                                            "Gene_size",
                                            "Percent_A",
                                            "Percent_T",
                                            "Percent_G",
                                            "Percent_C",
                                            "Nucleotide_Gene_Pos",
                                            "Nucleotide_chromosome_position",
                                            "Chromosome",
                                            "Reference_Nucleotide",
                                            "Alternative_Nucleotide",
                                            "Reference_Count",
                                            "Alternative_Count",
                                            "Tumor_origin",
                                            "Unique_patient_identifier",
                                            "Nucleotide_change_tally",
                                            "Nucleotide_mutation_rate",
                                            "Amino_acid_position",
                                            "Amino_acid_codon",
                                            "Codon_position",
                                            "Amino_acid_reference",
                                            "Amino_acid_alternative",
                                            "Amino_acid_mutation_rate",
                                            "Amino_acid_change_tally",
                                            "Gamma",
                                            "MAF_location",
                                            "Nucleotide_trinuc_context",
                                            "gene_level_synonymous_mutation_rate",
                                            "strand",
                                            "trinucs")

        mutation.data.to.add$Gene <- this.gene
        mutation.data.to.add$Gene_isoform <- RefCDS[[which(gene_list==this.gene)]]$protein_ID
        mutation.data.to.add$Gene_size <- length(myseq.split)
        mutation.data.to.add$Percent_A <- length(which(myseq.split=="A"))/length(myseq.split)
        mutation.data.to.add$Percent_T <- length(which(myseq.split=="T"))/length(myseq.split)
        mutation.data.to.add$Percent_G <- length(which(myseq.split=="G"))/length(myseq.split)
        mutation.data.to.add$Percent_C <- length(which(myseq.split=="C"))/length(myseq.split)
        mutation.data.to.add$Chromosome <- rep(this.gene_MAF[1,chr_column],length(mutation.data.to.add$Chromosome))
        mutation.data.to.add$gene_level_synonymous_mutation_rate <- gene_level_synonymous_mutation_rate
        mutation.data.to.add$strand <- this.strand

        mismatched.with.ref.genome <- NULL

        ####Tallying up mutations ----
        ##Loop that tallys up mutations in this gene, and their positions / affect on amino acid sequence.
        for(i in 1:nrow(this.gene_MAF)){

          #get the position within the gene and the amino acid code

          this.pos <- gene.pos[which(chrom.gene.pos==this.gene_MAF[i,pos_column])]

          this.AA.pos <- codon.numbers[this.pos]





          mutation.data.to.add$MAF_location[i] <- rownames(this.gene_MAF)[i]
          mutation.data.to.add$Nucleotide_chromosome_position[i] <- this.gene_MAF[i,pos_column]



          #what was the reference and the mutation
          if(this.strand=="-"){
            this.ref <- strand.switch[this.gene_MAF[i,ref_column],]
            this.mut <- strand.switch[this.gene_MAF[i,alt_column],]
          }else{
            this.ref <- this.gene_MAF[i,ref_column]
            this.mut <- this.gene_MAF[i,alt_column]
          }

          mutation.data.to.add$Reference_Nucleotide[i] <- this.ref
          mutation.data.to.add$Alternative_Nucleotide[i] <- this.mut
          mutation.data.to.add$Tumor_origin[i] <- this.gene_MAF[i,sample_ID_column]
          mutation.data.to.add$Unique_patient_identifier[i] <- this.gene_MAF[i,sample_ID_column]
          if(length(which(colnames(this.gene_MAF)=="t_ref_count"))>0){mutation.data.to.add$Reference_Count[i] <- this.gene_MAF$t_ref_count[i]}
          if(length(which(colnames(this.gene_MAF)=="t_alt_count"))>0){mutation.data.to.add$Alternative_Count[i] <- this.gene_MAF$t_alt_count[i]}

          #if this position exists (else, it wasn't in the isoform and tell us)...
          if(length(this.pos)>0){
            mutation.data.to.add$Nucleotide_Gene_Pos[i] <- this.pos
            mutation.data.to.add$Amino_acid_position[i] <- this.AA.pos
            mutation.data.to.add$Nucleotide_mutation_rate[i] <- normalized.mut.matrix[this.mut,this.pos]


            mutation.data.to.add$Nucleotide_trinuc_context[i] <- trinuc.pos[which(chrom.gene.pos==this.gene_MAF[i,pos_column])]

            #get which position within the codon the mutation occurred in
            if(this.pos%%3==1){pos <- 1};if(this.pos%%3==2){pos <- 2};if(this.pos%%3==0){pos <- 3}
            mutation.data.to.add$Codon_position[i] <- pos
            this.codon <- myseq.codons[this.AA.pos]
            mutation.data.to.add$Amino_acid_codon[i] <- this.codon
            this.codon.mutated <- unlist(strsplit(this.codon,split=""))
            this.codon.mutated[pos] <- this.mut
            this.codon.mutated <- paste(this.codon.mutated,collapse = "")


            this.codon.expanded <- unlist(strsplit(this.codon,split=""))

            this.AA.change <- AA_translations[which(AA_translations[,"Nucs"]==this.codon.mutated),"AA_short"]


            ##tally up the mutation in the tally matrices
            Nuc.tally.matrix[this.mut,this.pos] <- Nuc.tally.matrix[this.mut,this.pos]+1
            AA.tally.matrix[this.AA.change,this.AA.pos] <- AA.tally.matrix[this.AA.change,this.AA.pos]+1

            mutation.data.to.add$Amino_acid_reference[i] <- AA_translations[which(AA_translations[,"Nucs"]==this.codon)[1],"AA_letter"]
            mutation.data.to.add$Amino_acid_alternative[i] <- AA_translations[which(AA_translations[,"AA_short"]==this.AA.change)[1],"AA_letter"]

          }else{
            # print(paste("A mutation is outside the isoform of choice! The gene-specific MAF file position is: ",i,sep=""))

            if(length(which(outside.nucs$nuc_pos==this.gene_MAF[i,pos_column] & outside.nucs$change==this.mut))==0){
              open.spot <- which(is.na(outside.nucs$nuc_pos))[1]
              outside.nucs$nuc_pos[open.spot] <- this.gene_MAF[i,pos_column]
              outside.nucs$ref[open.spot] <- this.ref
              outside.nucs$change[open.spot] <- this.mut
              outside.nucs$tally[open.spot] <- 1

              #calculate the mutation rate at this sepecific site...
              if(this.strand=="-"){
                outside.sequence  <- rev(strsplit(as.character(BSgenome::getSeq(custom_BSgenome_selection,chromosome,start=this.gene_MAF[i,pos_column]-1,end=this.gene_MAF[i,pos_column]+1,strand="-")),split="")[[1]])
              }else{
                outside.sequence  <- strsplit(as.character(BSgenome::getSeq(custom_BSgenome_selection,chromosome,start=this.gene_MAF[i,pos_column]-1,end=this.gene_MAF[i,pos_column]+1,strand="+")),split="")[[1]]
              }

              mutation.data.to.add$Nucleotide_trinuc_context[i] <- paste(outside.sequence,collapse = "")

              #calculating the mutation rate for this gene
              if(outside.nucs$ref[open.spot]== "C" | outside.nucs$ref[open.spot]== "T"){
                if(outside.sequence[2]==outside.nucs$ref[open.spot]){
                  outside.nucs$mutation_rate[open.spot] <- (trinuc.mutation_data.num.this.gene[paste0(outside.sequence[2],outside.sequence[1],outside.nucs$change[open.spot],outside.sequence[3],collapse = ""),"proportion"]/mean.nuc.rate.for.this.gene)*gene_level_synonymous_mutation_rate
                  mutation.data.to.add$trinucs[i] <- paste0(outside.sequence[2],outside.sequence[1],outside.nucs$change[open.spot],outside.sequence[3],collapse = "")
                }else{
                  #This means that the reference nucleotide in the MAF does not match the reference genome position.
                  message(paste("The reference allele given in the MAF file does not match the reference genome!\nThis gene: ",this.gene," \nThis position in the gene-specific MAF file: ",i)) #Let us know
                  mismatched.with.ref.genome <- c(mismatched.with.ref.genome,i)
                }
              }else{
                if(outside.sequence[2]==outside.nucs$ref[open.spot]){
                  outside.nucs$mutation_rate[open.spot] <- (trinuc.mutation_data.num.this.gene[paste0(flip.function(outside.sequence[2]),flip.function(outside.sequence[3]),flip.function(outside.nucs$change[open.spot]),flip.function(outside.sequence[1]),collapse = ""),"proportion"]/mean.nuc.rate.for.this.gene)*gene_level_synonymous_mutation_rate
                  mutation.data.to.add$trinucs[i] <- paste0(flip.function(outside.sequence[2]),flip.function(outside.sequence[3]),flip.function(outside.nucs$change[open.spot]),flip.function(outside.sequence[1]),collapse = "")
                }else{
                  #This means that the reference nucleotide in the MAF does not match the reference genome position.
                  message(paste("The reference allele given in the MAF file does not match the reference genome!\nThis gene: ",this.gene," \nThis position in the gene-specific MAF file: ",i)) #Let us know
                  mismatched.with.ref.genome <- c(mismatched.with.ref.genome,i)
                }
              }
              mutation.data.to.add$Nucleotide_mutation_rate[i] <- outside.nucs$mutation_rate[open.spot]
            }else{
              outside.nucs$tally[which(outside.nucs$nuc_pos==this.gene_MAF[i,pos_column] & outside.nucs$change==this.mut)] <- outside.nucs$tally[which(outside.nucs$nuc_pos==this.gene_MAF[i,pos_column] & outside.nucs$change==this.mut)]+1

            }

            # mutation.data.to.add$Nucleotide_mutation_rate[i] <-     #need to calculate for this trinucleotide

          }
        }
        if(length(which(is.na(outside.nucs$nuc_pos)))>0){
          outside.nucs <- outside.nucs[-which(is.na(outside.nucs$nuc_pos)),] #reducing data frame size
        }
        if(length(mismatched.with.ref.genome)>0){
          mutation.data.to.add <- mutation.data.to.add[-mismatched.with.ref.genome,]
        }





        #codon.numbers

        mut.spots.nuc <- which(Nuc.tally.matrix!=0,arr.ind = T)
        if(nrow(mut.spots.nuc)==0 & nrow(outside.nucs)==0){
          print("No mutations to tally!")
          next
        }

        if(nrow(mut.spots.nuc)>0){
          for(i in 1:nrow(mut.spots.nuc)){
            these.pos <- which(mutation.data.to.add$Nucleotide_Gene_Pos==mut.spots.nuc[i,2] & mutation.data.to.add$Alternative_Nucleotide==rownames(mut.spots.nuc)[i])
            mutation.data.to.add$Nucleotide_change_tally[these.pos] <- length(these.pos)
          }
        }


        #TODO, do not need to generate entire AA.gamma.matrix, may be wasting memory and computation time.
        mut.spots <- which(AA.tally.matrix!=0,arr.ind=T)
        if(nrow(mut.spots)==0 & nrow(outside.nucs)==0){
          print("No mutations to tally!")
          next
        }



        #tallying up the changes!
        if(nrow(mut.spots)>0){
          for(i in 1:nrow(mut.spots)){
            these.pos <- which(mutation.data.to.add$Amino_acid_position==mut.spots[i,2] & mutation.data.to.add$Amino_acid_alternative==AA_translations[which(AA_translations[,"AA_short"]==rownames(mut.spots)[i])[1],"AA_letter"])

            mutation.data.to.add$Amino_acid_change_tally[these.pos] <- length(these.pos)

            mutation.data.to.add$Amino_acid_mutation_rate[these.pos] <- AA.mut.matrix[mut.spots[i,1],mut.spots[i,2]]

            #TODO: had to remove this for size constraint reasons.

            # mutation.data.to.add$trinucs[these.pos] <- all_genes$trinuc.comp[[which(all_gene_trinuc_gene_list==this.gene)]]$trinuc_matrix[mut.spots[i,1],mut.spots[i,2]]


            # mutation.data.to.add$trinucs[these.pos] <- AA.matrix.trinuc[mut.spots[i,1],mut.spots[i,2]]




            if(AA.mut.matrix[mut.spots[i,1],mut.spots[i,2]]==0){message(paste("The mutation rate for this mutation is listed as ZERO"));print(mut.spots[i,])
            }else{
              AA.gamma.matrix[mut.spots[i,1],mut.spots[i,2]] <- lambda.calc(n.one = AA.tally.matrix[mut.spots[i,1],mut.spots[i,2]],n.zero = (tumor.number-AA.tally.matrix[mut.spots[i,1],mut.spots[i,2]]))/AA.mut.matrix[mut.spots[i,1],mut.spots[i,2]]
              # newton.method(p = AA.tally.matrix[mut.spots[i,1],mut.spots[i,2]],n = tumor.number,mutation_rate = AA.mut.matrix[mut.spots[i,1],mut.spots[i,2]])
            }

            mutation.data.to.add$Gamma[these.pos] <-  AA.gamma.matrix[mut.spots[i,1],mut.spots[i,2]]
          }
        }
        if(nrow(outside.nucs)>0){
          outside.nucs$selection_intensity <- NA
          outside.nucs$selection_no_epistasis <- NA
          for(i in 1:nrow(outside.nucs)){

            outside.nucs$selection_intensity[i] <- lambda.calc(n.one = outside.nucs$tally[i],n.zero = (tumor.number - length(unique(this.gene_MAF[,sample_ID_column]))))/outside.nucs$mutation_rate[i]
            outside.nucs$selection_no_epistasis[i] <- lambda.calc(n.one = outside.nucs$tally[i],n.zero = (tumor.number - outside.nucs$tally[i]))/outside.nucs$mutation_rate[i]

            mutation.data.to.add$Nucleotide_change_tally[which(mutation.data.to.add$Nucleotide_chromosome_position==outside.nucs$nuc_pos[i])] <- outside.nucs$tally[i]
          }

        }

        #Add the gamma assuming complete epistasis calculation to the complete_mutation_data dataframe
        mutation.data.to.add$selection_intensity <- NA

        for(i in 1:nrow(mutation.data.to.add)){
          if(is.na(mutation.data.to.add$Amino_acid_position[i])){
            mutation.data.to.add$selection_intensity[i] <- lambda.calc(n.one = mutation.data.to.add$Nucleotide_change_tally[i],n.zero = (tumor.number - length(unique(this.gene_MAF[,sample_ID_column]))))/mutation.data.to.add$Nucleotide_mutation_rate[i]
          }else{
            mutation.data.to.add$selection_intensity[i] <- lambda.calc(n.one = mutation.data.to.add$Amino_acid_change_tally[i],n.zero = (tumor.number - length(unique(this.gene_MAF[,sample_ID_column]))))/mutation.data.to.add$Amino_acid_mutation_rate[i]
          }
        }


        this.gene.gamma <- which(!is.na(AA.gamma.matrix),arr.ind=T)
        added.rows <- as.data.frame(matrix(data=NA,nrow=(nrow(this.gene.gamma)+nrow(outside.nucs)),ncol=ncol(all.muts)),header=T)
        colnames(added.rows) <- colnames(all.muts)
        added.rows$Gene <- this.gene
        added.rows$gene_AA_size <- length(myseq.codons)
        added.rows$dndscv_p <- dndscv_siggenes$p[which(dndscv_siggenes$gene==this.gene)]
        added.rows$dndscv_q <- dndscv_siggenes$q[which(dndscv_siggenes$gene==this.gene)]
        added.rows$Prop_of_tumors_with_this_gene_mutated <- length(unique(this.gene_MAF[,sample_ID_column]))/tumor.number

        counter <- 1
        if(nrow(this.gene.gamma)>0){
          for(e in 1:nrow(this.gene.gamma)){
            added.rows$AA_Pos[e] <- this.gene.gamma[e,2]
            added.rows$AA_Ref[e] <- AA.code[added.rows$AA_Pos[e]]
            added.rows$AA_Change[e] <- AA_translations[which(AA_translations[,"AA_short"]==rownames(AA.gamma.matrix)[this.gene.gamma[e,1]])[1],"AA_letter"]
            added.rows$gamma[e] <- AA.gamma.matrix[this.gene.gamma[e,1],this.gene.gamma[e,2]]
            added.rows$mutation_rate[e] <- AA.mut.matrix[this.gene.gamma[e,1],this.gene.gamma[e,2]]
            added.rows$prevalence_among_tumors[e] <- AA.tally.matrix[this.gene.gamma[e,1],this.gene.gamma[e,2]]
            added.rows$Prop_tumors_with_specific_mut[e] <- AA.tally.matrix[this.gene.gamma[e,1],this.gene.gamma[e,2]]/tumor.number
            added.rows$selection_intensity[e] <- lambda.calc(n.one = added.rows$prevalence_among_tumors[e],n.zero = (tumor.number - length(unique(this.gene_MAF[,sample_ID_column]))))/added.rows$mutation_rate[e]
            counter <- counter + 1
          }
        }
        if(nrow(outside.nucs)>0){
          for(i in 1:nrow(outside.nucs)){
            added.rows[counter,c("Nucleotide_position","Nuc_Ref","Nuc_Change","prevalence_among_tumors","mutation_rate","selection_intensity","gamma")] <-
              outside.nucs[i,c("nuc_pos","ref","change","tally","mutation_rate","selection_intensity","selection_no_epistasis")]
            added.rows[counter,"Prop_tumors_with_specific_mut"] <- added.rows[counter,"prevalence_among_tumors"]/tumor.number
            counter <- counter + 1
          }
        }
        # added.rows <- complete.epi.function(mutationdf = added.rows,samplenumber = tumor.number)

        # print(added.rows)
        all.muts <- rbind(all.muts,added.rows)

        mutation.data <- rbind(mutation.data,mutation.data.to.add)
      }else{
        #This is what happens if you already have the mutation data (isoform, trinucleotide distribution per gene, trinucleotide context per substitution)
        #output_from_mainMAF should have all possible mutations in the bootstrap MAFs

        if(length(which(output_from_mainMAF$complete_mutation_data$Gene==this.gene))==0){
          message(paste("There is no previous selection data for ",this.gene))
          next
        }


        main_MAF.this.gene <- output_from_mainMAF$complete_mutation_data[which(output_from_mainMAF$complete_mutation_data$Gene==this.gene),]
        this.strand <- main_MAF.this.gene$strand[1]

        trinuc.mutation_data.num.this.gene <- trinuc.mutation_data.num
        # trinuc.data.this.gene <- trinuc.mutation_data
        # trinuc.data.this.gene$total_count <- 0
        trinuc.mutation_data.num.this.gene[,"total_count"] <- output_from_mainMAF$trinuc_counts[,this.gene]

        to.mean <-NULL
        for(w in 1:32){
          # to.mean[w] <- trinuc.data.this.gene$total_count[w*3] * sum(trinuc.data.this.gene$proportion[c(w*3,(w*3)-1,(w*3)-2)])
          to.mean[w] <- trinuc.mutation_data.num.this.gene[w*3,"total_count"] * sum(trinuc.mutation_data.num.this.gene[c(w*3,(w*3)-1,(w*3)-2),"proportion"])
        }
        mean.nuc.rate.for.this.gene <- (sum(to.mean)/(sum(trinuc.mutation_data.num.this.gene[,"total_count"])/3))

        # sometimes the same gene name is on different chromosomes. This filters to the isoform of choice in the previous analysis
        this.gene_MAF <- this.gene_MAF[which(this.gene_MAF[,chr_column]==main_MAF.this.gene[1,chr_column]),]

        # Need to cut out substitutions that are inconsistent between the MAF file and the reference genome. These were recorded in the previous run.
        if(length(which(output_from_mainMAF$cut.mutations[,gene_ID_column]==this.gene))>0){
          this.gene_MAF <- this.gene_MAF[-which(this.gene_MAF[,gene_ID_column] %in% output_from_mainMAF$cut.mutations[,gene_ID_column] & this.gene_MAF[,pos_column] %in% output_from_mainMAF$cut.mutations[,pos_column] & this.gene_MAF[,chr_column] %in% output_from_mainMAF$cut.mutations[,chr_column]),]
        }


        if(nrow(this.gene_MAF)==0){
          print(paste("All the data for this gene is from another chromosome, after cleaning based on the isoform picked for the original analysis"))
          next
        }

        #If it is on the - strand, we need to flip the variants to that strand so it is consistent with previous output. MAFs always contain "+" strand data
        if(this.strand=="-"){
          this.gene_MAF[,ref_column] <- unlist(lapply(X = this.gene_MAF[,ref_column],FUN = flip.function))
          this.gene_MAF[,alt_column] <- unlist(lapply(X = this.gene_MAF[,alt_column],FUN = flip.function))
        }



        gene_level_synonymous_mutation_rate <- mut_rates[this.gene]



        added.rows <- as.data.frame(matrix(data=NA,nrow=nrow(this.gene_MAF),ncol=ncol(all.muts)),header=T)
        colnames(added.rows) <- colnames(all.muts)
        added.rows$Gene <- this.gene
        added.rows$dndscv_p <- dndscv_siggenes$p[which(dndscv_siggenes$gene==this.gene)]
        added.rows$dndscv_q <- dndscv_siggenes$q[which(dndscv_siggenes$gene==this.gene)]
        added.rows$Prop_of_tumors_with_this_gene_mutated <- length(unique(this.gene_MAF[,sample_ID_column]))/tumor.number

        for(mafRow in 1:nrow(this.gene_MAF)){


          if(!is.na(main_MAF.this.gene$Amino_acid_position[which(main_MAF.this.gene$Nucleotide_chromosome_position==this.gene_MAF[mafRow,pos_column])[1]])){ #within the exon
            if(length(which(added.rows$AA_Pos == main_MAF.this.gene$Amino_acid_position[which(main_MAF.this.gene$Nucleotide_chromosome_position==this.gene_MAF[mafRow,pos_column])[1]] & added.rows$AA_Change == main_MAF.this.gene$Amino_acid_alternative[which(main_MAF.this.gene$Nucleotide_chromosome_position==this.gene_MAF[mafRow,pos_column] & main_MAF.this.gene$Alternative_Nucleotide==this.gene_MAF[mafRow,alt_column])[1]]))==0){ #If this mutation isn't already stored in the output
              this.row <- which(is.na(added.rows$mutation_rate))[1]

              added.rows$AA_Pos[this.row] <- main_MAF.this.gene$Amino_acid_position[which(main_MAF.this.gene$Nucleotide_chromosome_position==this.gene_MAF[mafRow,pos_column])[1]]
              added.rows$AA_Ref[this.row] <- main_MAF.this.gene$Amino_acid_reference[which(main_MAF.this.gene$Nucleotide_chromosome_position==this.gene_MAF[mafRow,pos_column])[1]]
              this.AA.change <- main_MAF.this.gene$Amino_acid_alternative[which(main_MAF.this.gene$Nucleotide_chromosome_position==this.gene_MAF[mafRow,pos_column] & main_MAF.this.gene$Alternative_Nucleotide == this.gene_MAF[mafRow,alt_column])[1]]
              #Need the rows in the main_MAF.this.gene that are also this AA pos and change, and their unique nucleotide positions and nucleotide changes, and then
              #if the this.gene_MAF file has these as well, to tally up!
              added.rows$AA_Change[this.row] <- this.AA.change


              nuc.changes <- unlist(strsplit(main_MAF.this.gene$trinucs[which(main_MAF.this.gene$Nucleotide_chromosome_position==this.gene_MAF[mafRow,pos_column] & main_MAF.this.gene$Alternative_Nucleotide == this.gene_MAF[mafRow,alt_column])[1]],split=":"))




              added.rows$mutation_rate[this.row] <- sum(trinuc.mutation_data.num.this.gene[nuc.changes,"proportion"])/mean.nuc.rate.for.this.gene*gene_level_synonymous_mutation_rate


              #now, need to find all spots in the this.geneMAF that correspond to this AA change.
              #Need to find all mutation types in the main_MAF file that correspond to this AA change, and count them in the this.geneMAF file
              this.main_MAF.pos <- which(main_MAF.this.gene$Nucleotide_chromosome_position==this.gene_MAF[mafRow,pos_column] & main_MAF.this.gene$Alternative_Nucleotide == this.gene_MAF[mafRow,alt_column])[1]
              all.AA.mut.positions <- which(main_MAF.this.gene$Amino_acid_alternative==main_MAF.this.gene$Amino_acid_alternative[this.main_MAF.pos] & main_MAF.this.gene$Amino_acid_position==main_MAF.this.gene$Amino_acid_position[this.main_MAF.pos])

              #all nuc positions and changes that resulted in this AA pos and AA change
              these.nuc.pos <- main_MAF.this.gene$Nucleotide_chromosome_position[all.AA.mut.positions]
              these.nuc.change <- main_MAF.this.gene$Alternative_Nucleotide[all.AA.mut.positions]

              this.gene_MAF.nucs <- which(this.gene_MAF[,pos_column] %in% these.nuc.pos)
              this.gene_MAF.change <- which(this.gene_MAF[,alt_column] %in% these.nuc.change)

              added.rows$prevalence_among_tumors[this.row] <- length(intersect(this.gene_MAF.nucs,this.gene_MAF.change))

              added.rows$selection_intensity[this.row] <- lambda.calc(n.one = added.rows$prevalence_among_tumors[this.row],n.zero = (tumor.number-length(unique(this.gene_MAF[,sample_ID_column]))))/added.rows$mutation_rate[this.row]
              added.rows$gamma[this.row] <- lambda.calc(n.one = added.rows$prevalence_among_tumors[this.row],n.zero = (tumor.number-added.rows$prevalence_among_tumors[this.row]))/added.rows$mutation_rate[this.row]

            } #if it is, just skip it. Should have stored all the information on previous loops



          }else{ #outside the exon

            #only do the calculation if it is not already present in the added.rows dataframe
            if(length(which(added.rows$Nucleotide_position==this.gene_MAF[mafRow,pos_column] & added.rows$Nuc_Change==this.gene_MAF[mafRow,alt_column]))==0){

              this.row <- which(is.na(added.rows$mutation_rate))[1]
              added.rows$Nucleotide_position[this.row] <- this.gene_MAF[mafRow,pos_column]
              added.rows$Nuc_Ref[this.row] <- this.gene_MAF[mafRow,ref_column]
              added.rows$Nuc_Change[this.row] <- this.gene_MAF[mafRow,alt_column]

              # added.rows$mutation_rate[this.row] <- (main_MAF.this.gene$Nucleotide_mutation_rate[which(main_MAF.this.gene$Nucleotide_chromosome_position==added.rows$Nucleotide_position[this.row] & main_MAF.this.gene$Alternative_Nucleotide==added.rows$Nuc_Change[this.row])[1]]*gene_level_synonymous_mutation_rate)/(main_MAF.this.gene$gene_level_synonymous_mutation_rate[1])
              nuc.changes <- unlist(strsplit(main_MAF.this.gene$trinucs[which(main_MAF.this.gene$Nucleotide_chromosome_position==this.gene_MAF[mafRow,pos_column] & main_MAF.this.gene$Alternative_Nucleotide == this.gene_MAF[mafRow,alt_column])[1]],split=":"))
              added.rows$mutation_rate[this.row] <- sum(trinuc.mutation_data.num.this.gene[nuc.changes,"proportion"])/mean.nuc.rate.for.this.gene*gene_level_synonymous_mutation_rate

              added.rows$prevalence_among_tumors[this.row] <- length(which(this.gene_MAF[,pos_column]==added.rows$Nucleotide_position[this.row] & this.gene_MAF[,alt_column]==added.rows$Nuc_Change[this.row]))

              added.rows$selection_intensity[this.row] <- lambda.calc(n.one = added.rows$prevalence_among_tumors[this.row],n.zero = (tumor.number-length(unique(this.gene_MAF[,sample_ID_column]))))/added.rows$mutation_rate[this.row]

              added.rows$gamma[this.row] <- lambda.calc(n.one = added.rows$prevalence_among_tumors[this.row],n.zero = (tumor.number-added.rows$prevalence_among_tumors[this.row]))/added.rows$mutation_rate[this.row]

            }
          }



        }
        added.rows <- added.rows[which(!is.na(added.rows$mutation_rate)),]
        added.rows$Prop_tumors_with_specific_mut <- added.rows$prevalence_among_tumors/tumor.number
        all.muts <- rbind(all.muts,added.rows)
      }
    }else{
      print(paste("The gene ",this.gene," is not in the MAF file!",sep=""))
    }

    if(zzz%%1000==0){message(paste("Gene number: ",zzz," out of ",length(genes_for_analysis)))} #If zzz is a multiple of 1000 print the progress

  }
  if(exists('nuc.mutation.rate.vec')){ #this catch means the function generated the gene sequence. If not, just store what we have
    output.list <- list(all_mutations=all.muts,
                        complete_mutation_data=mutation.data,
                        nucleotide_mutation_vec=nuc.mutation.rate.vec,
                        nucleotide_tally=Nuc.tally.matrix,
                        amino_acid_tally=AA.tally.matrix,
                        amino_acid_mutation_rates=AA.mut.matrix,
                        nucleotide_mutation_rates=normalized.mut.matrix,
                        myseqsplit=myseq.split,
                        # isoform.list=isoform.list,
                        trinuc_counts=trinuc.count.matrix,
                        cut.mutations=mutations.to.cut[-is.na(mutations.to.cut[,gene_ID_column]),])
  }else{

    all.muts$identifier <- NA
    for(i in 1:nrow(all.muts)){
      if(!(is.na(all.muts$AA_Pos[i]))){
        all.muts$identifier[i] <- paste0(all.muts$Gene[i]," ",all.muts$AA_Ref[i],all.muts$AA_Pos[i],all.muts$AA_Change[i])
      }else{
        all.muts$identifier[i] <- paste0(all.muts$Gene[i]," ",all.muts$Nuc_Ref[i],all.muts$Nucleotide_position[i],all.muts$Nuc_Change[i])
      }
    }

    output.list <- list(all_mutations=all.muts,
                        complete_mutation_data=mutation.data)

  }


  return(output.list)

}

