#' buildref
#' 
#' Function to build a RefCDS object from a reference genome and a table of transcripts. The RefCDS object has to be precomputed for any new species or assembly prior to running dndscv. This function generates an .rda file that needs to be input into dndscv using the refdb argument. Note that when multiple CDS share the same gene name (second column of cdsfile), the longest coding CDS will be chosen for the gene. CDS with ambiguous bases (N) will not be considered.
#' @importFrom crayon silver
#' @details Based on the buildref function in Inigo Martincorena's package dNdScv.
#' @param cds_data data from biomaRt for the genome assembly of interest
#' @param genome a BSgenome object corresponding to the genome assembly of interest
#' @param cores how many cores to use for parallel computations (requires parallel package)
#' @param numcode NCBI genetic code number (default = 1; standard genetic code). To see the list of genetic codes supported use: ? seqinr::translate
#' @param excludechrs Vector or string with chromosome names to be excluded from the RefCDS object (default: no chromosome will be excluded). The mitochondrial chromosome should be excluded as it has different genetic code and mutation rates, either using the excludechrs argument or not including mitochondrial transcripts in cdsfile.
#' @param onlychrs Vector of valid chromosome names (default: all chromosomes will be included)
#' @param useids Combine gene IDs and gene names (columns 1 and 2 of the input table) as long gene names (default = F)
#' 
#' @export

build_RefCDS = function(cds_data, genome, numcode = 1, excludechrs = NULL, onlychrs = NULL, useids = F, cores=1) {
  message("[Step 1/3] Identifying complete transcripts...")
  # Valid chromosomes and reference CDS per gene
  reftable = cds_data
  GenomeInfoDb::seqlevelsStyle(genome) = "NCBI" # assuming this doesn't crash on any genomes
  
  # Support for alternate genetic codes not yet implemented
  if(numcode != 1) {
    stop("Currently only the standard genetic code is supported.")
  }
  
  # change column names
  colnames(reftable) = c("gene.id","cds.id","chr","strand","gene.name","length","chr.coding.start","chr.coding.end","cds.start","cds.end")
  reftable[,6:10] = suppressWarnings(lapply(reftable[,6:10], as.numeric)) # Convert to numeric
  
  # Checking for systematic absence of gene names (it happens in some BioMart inputs)
  longname = paste(reftable$gene.id, reftable$gene.name, sep=":") # Gene name combining the Gene stable ID and the Associated gene name
  if (useids==T) {
    reftable$gene.name = longname # Replacing gene names by a concatenation of gene ID and gene name
  }
  if (length(unique(reftable$gene.name))<length(unique(longname))) {
    warning(sprintf("%0.0f unique gene IDs (column 1) found. %0.0f unique gene names (column 2) found. Consider combining gene names and gene IDs or replacing gene names by gene IDs to avoid losing genes (see useids argument in ? buildref)",length(unique(reftable$gene.id)),length(unique(reftable$gene.name))))
  }
  
  # Reading chromosome names in the fasta file using its index file if it already exists or creating it if it does not exist. The index file is also used by scanFa later.
  

  validchrs = as.character(GenomicRanges::seqnames(genome))
  validchrs = setdiff(validchrs, excludechrs)
  if (length(onlychrs)>0) {
    validchrs = validchrs[validchrs %in% onlychrs]
  }
  
  # Restricting to chromosomes present in both the genome file and the CDS table
  if (any(validchrs %in% unique(reftable$chr))) {
    validchrs = validchrs[validchrs %in% unique(reftable$chr)]
  }
  
  
  # Removing genes that fall partially or completely outside of the available chromosomes/contigs
  reftable = reftable[reftable[,1]!="" & reftable[,2]!="" & reftable[,5]!="" & !is.na(reftable[,7]) & !is.na(reftable[,8]),] # Removing invalid entries
  reftable = reftable[which(as.character(reftable$chr) %in% validchrs),] # Only valid chromosomes
  
  transc_gr = GenomicRanges::GRanges(reftable$chr, IRanges::IRanges(reftable$chr.coding.start,reftable$chr.coding.end))
  chrs_gr = GenomicRanges::GRanges(GenomicRanges::seqinfo(genome))
  ol = as.data.frame(GenomicRanges::findOverlaps(transc_gr, chrs_gr, type="within", select="all"))
  
  # Issuing an error if any transcript falls outside of the limits of a chromosome. Possibly due to a mismatch between the assembly used for the reference table and the reference genome.
  if (length(unique(ol[,1])) < nrow(reftable)) {
    stop(sprintf("Aborting buildref. %0.0f rows in cdsfile have coordinates that fall outside of the corresponding chromosome length. Please ensure that you are using the same assembly for the cdsfile and genomefile",nrow(reftable)-length(unique(ol[,1]))))
  }
  
  reftable = reftable[unique(ol[,1]),] 
  
  # Identifying genes starting or ending at the ends of a chromosome/contig
  # Because buildref and dndscv need to access the base before and after each coding position, genes overlapping the ends
  # of a contig will be trimmed by three bases and a warning will be issued listing those genes.
  
  fullcds = intersect(reftable$cds.id[reftable$cds.start==1], reftable$cds.id[reftable$cds.end==reftable$length]) # List of complete CDS
  
  ol_start = as.data.frame(GenomicRanges::findOverlaps(transc_gr, chrs_gr, type="start", select="all"))[,1] # Genes overlapping contig starts
  if (any(ol_start)) {
    reftable[ol_start,"chr.coding.start"] = reftable[ol_start,"chr.coding.start"] + 3 # Truncate the first 3 bases
    reftable[ol_start,"cds.start"] = reftable[ol_start,"cds.start"] + 3 # Truncate the first 3 bases
  }
  
  ol_end = as.data.frame(GenomicRanges::findOverlaps(transc_gr, chrs_gr, type="end", select="all"))[,1] # Genes overlapping contig starts
  if (any(ol_end)) {
    reftable[ol_end,"chr.coding.end"] = reftable[ol_end,"chr.coding.end"] - 3 # Truncate the first 3 bases
    reftable[ol_end,"cds.end"] = reftable[ol_end,"cds.end"] - 3 # Truncate the first 3 bases
  }
  
  if (any(c(ol_start,ol_end))) {
    warning(sprintf("The following genes were found to start or end at the first or last base of their contig. Since dndscv needs trinucleotide contexts for all coding bases, codons overlapping ends of contigs have been trimmed. Affected genes: %s.", paste(reftable[unique(c(ol_start,ol_end)),"gene.name"], collapse=", ")))
  }
  
  
  # Selecting the longest complete CDS for every gene (required when there are multiple alternative transcripts per unique gene name)
  reftable = reftable[reftable$cds.id %in% fullcds, ] # Filter to just complete transcripts (some will have missing exons)
  
  # Select on reftable row per CDS to create CDS "dictionary"
  cds_table = unique(reftable[,c(1,2,5,6)])
  cds_table = cds_table[order(cds_table$gene.name, -cds_table$length), ] # Sorting CDS from longest to shortest
  cds_table = cds_table[(cds_table$length %% 3)==0, ] # Removing CDS of length not multiple of 3
  gene_list = unique(cds_table$gene.name)
  
  reftable = reftable[order(reftable$chr, reftable$chr.coding.start), ]

  cds_split = split(reftable, f=reftable$cds.id)
  gene_split = split(cds_table, f=cds_table$gene.name)
  
  message("[Step 2/3] Examining exons and splice sites of complete transcripts...")
  remaining_genes = gene_list
  which_transcript = 1
  num_transcript_by_gene = sapply(gene_split, nrow)
  invalid_genes = character()
  ref_env = new.env()
  while(length(remaining_genes) > 0) {
    message(silver(paste0("\tPass ", which_transcript, ":")))
    remaining_gene_info = gene_split[remaining_genes]
    transcript_ids = sapply(remaining_gene_info, function(x) x[which_transcript,"cds.id"])
    cds_split_subset = cds_split[transcript_ids]
    reftable_subset = reftable[reftable$cds.id %in% transcript_ids,]
    reftable_subset$strand[reftable_subset$strand == "-1"] = "-"
    reftable_subset$strand[reftable_subset$strand == "1"] = "+"
    
    grs = GenomicRanges::makeGRangesFromDataFrame(reftable_subset, keep.extra.columns = T, 
                                                  seqinfo = GenomicRanges::seqinfo(genome), seqnames.field = "chr", 
                                                  start.field = "chr.coding.start", end.field = "chr.coding.end", strand.field = "strand")
    GenomicRanges::start(grs) = GenomicRanges::start(grs) - 1
    GenomicRanges::end(grs) = GenomicRanges::end(grs) + 1
    
    # for efficiency, start with list object, convert to GRangesList for getSeq, then environment for get_spl_info
    grs_by_transcript = split(grs, f= grs$cds.id)
    
    # reverse granges exon order on negative-strand transcripts
    transcripts = names(grs_by_transcript)
    chrs = sapply(cds_split_subset[transcripts], function(x) x$chr[1])
    strands = sapply(cds_split_subset[transcripts], function(x) x$strand[1])
    strands[strands == "-1"] = "-"
    strands[strands == "1"] = "+"
    grs_by_transcript[strands == '-'] = lapply(grs_by_transcript[strands == '-'], rev)
    
    # Get sequences
    message(silver("\tGetting exonic DNA sequences..."))
    grangeslist = GenomicRanges::GRangesList(grs_by_transcript)
    
    padded_cds_seqs = list2env(as.list(BSgenome::getSeq(genome, grangeslist)))
    
    grs_by_transcript = list2env(grs_by_transcript)
    
    
    # Essential splice sites are 1,2 bp upstream (3') and 1,2,5 bp downstream (5')
    ## The inputs to this function were padded +/- 1 base, so we add/subtract 0,1 and 0,1,4 to get the same positions
    get_spl_info = function(pid) {
      gr = grs_by_transcript[[pid]]
      strand = strands[[pid]]
      chr = chrs[[pid]]
      # extend +/- strand positions accordingly
      if(strand == "+") {
        spl3prime = GenomicRanges::start(gr[-1]) # first exon has no 3' splice site
        spl5prime = GenomicRanges::end(gr[-length(gr)]) # last exon has no 5' splice site
        splpos = unique(sort(c(spl5prime, spl5prime+1, spl5prime+4, spl3prime, spl3prime-1)))
        
      } else {
        spl3prime = GenomicRanges::end(gr[-1]) # first exon has no 3' splice site
        spl5prime = GenomicRanges::start(gr[-length(gr)]) # last exon has no 5' splice site
        splpos = unique(sort(c(spl5prime, spl5prime-1, spl5prime-4, spl3prime, spl3prime+1)))
      }
      
      # empty GRanges object if there are no splice sites (happens for 1-exon transcripts)
      if (length(splpos) == 0) {
        grange = GenomicRanges::GRanges()
      } else {
        ranges = IRanges::IRanges(start = splpos-1, end = splpos + 1)
        if (strand == "-") {
          ranges = rev(ranges)
        }
        grange = GenomicRanges::GRanges(seqnames = chr, ranges = ranges, strand = strand)
      }
      return(list(splpos, grange))
    }
    
    transcript_names = names(grs_by_transcript)
    message(silver("\tCalculating splice site info..."))
    spl_info = pbapply::pblapply(transcript_names, get_spl_info, cl = cores)
    message(silver("\tGetting splice site DNA sequences..."))
    splpos_by_transcript = lapply(spl_info, function(x) x[[1]])
    splrange_by_transcript = GenomicRanges::GRangesList(lapply(spl_info, function(x) x[[2]]))
    names(splpos_by_transcript) = transcript_names
    names(splrange_by_transcript) = transcript_names
    padded_spl_seqs = list2env(as.list(BSgenome::getSeq(genome, splrange_by_transcript)))
    genes_for_next_time = character()
    
    
    for(gene in remaining_genes) {
      info = remaining_gene_info[[gene]]
      pid = info$cds.id[which_transcript]
      
      cds = cds_split_subset[[pid]]
      strand = cds[1,4]
      chr = cds[1,3]
      # get grange with extra base on both ends of all ranges
      cdsseq_padded = padded_cds_seqs[[pid]]
      padded_widths = Biostrings::width(cdsseq_padded)
      
      # exclude the 1bp padding at each end of each sequence, and paste all together to get full cds sequence
      cdsseq = unlist(Biostrings::subseq(cdsseq_padded, start = 2, end = padded_widths - 1))
      
      # similar idea to get sequences with upstream/downstream bases
      cdsseq1up = unlist(Biostrings::subseq(cdsseq_padded, start = 1, end = padded_widths - 2))
      cdsseq1down = unlist(Biostrings::subseq(cdsseq_padded, start = 3, end = padded_widths))
      
      
      # Get AA sequence (can add alternative translation codes later; looks like Biostrings supports them)
      pseq = as.character(Biostrings::translate(cdsseq))
      
      # A valid CDS has been found if this test passes
      # No stop codons ("*") in pseq excluding the last codon, and no "N" nucleotides in cdsseq
      if(! grepl("\\*.", pseq) && Biostrings::countPattern("N", cdsseq) == 0) {
        # Essential splice sites
        splpos = splpos_by_transcript[[pid]] # Essential splice sites
        if (length(splpos)>0) { # CDSs with a single exon do not have splice sites
          
          # Obtaining the splicing sequences and the coding and splicing sequence contexts
          # get splice ranges +/- 1 base
          splseq_padded = padded_spl_seqs[[pid]]
          padded_widths = Biostrings::width(splseq_padded)
          
          # exclude the 1bp padding at each end of each sequence, and paste all together to get full cds sequence
          splseq = unlist(Biostrings::subseq(splseq_padded, start = 2, end = padded_widths - 1))
          
          # similar idea to get sequences with upstream/downstream bases
          splseq1up = unlist(Biostrings::subseq(splseq_padded, start = 1, end = padded_widths - 2))
          splseq1down = unlist(Biostrings::subseq(splseq_padded, start = 3, end = padded_widths))
        }
        
        # Annotating the CDS in the RefCDS database
        gene_name = gene
        gene_id = info$gene.id[which_transcript]
        protein_id = pid
        CDS_length = info$length[which_transcript]
        chr = cds[1,3]
        strand = strand
        intervals_cds = unname(as.matrix(cds[,7:8]))
        intervals_splice = splpos
        seq_cds = cdsseq
        seq_cds1up = cdsseq1up
        seq_cds1down = cdsseq1down
        
        results = list(gene_name = gene_name, gene_id = gene_id, 
                       protein_id = protein_id, CDS_length = CDS_length, chr = chr, 
                       strand = strand, intervals_cds = intervals_cds,
                       intervals_splice = intervals_splice, seq_cds = seq_cds, 
                       seq_cds1up = seq_cds1up, seq_cds1down = seq_cds1down)
        
        if (length(splpos)>0) { # If there are splice sites in the gene
          seq_splice = splseq
          seq_splice1up = splseq1up
          seq_splice1down = splseq1down
          results = c(results, list(seq_splice = seq_splice, seq_splice1up = seq_splice1up, seq_splice1down = seq_splice1down))
        }
        ref_env[[gene_name]] = results
      } else {
        if(num_transcript_by_gene[[gene]] > which_transcript) {
          genes_for_next_time = c(genes_for_next_time, gene)
        } else {
          invalid_genes = c(invalid_genes, gene)
        }
      }
    }
    remaining_genes = genes_for_next_time
    which_transcript = which_transcript + 1
  }
  
  ## 3. L matrices: number of synonymous, missense, nonsense and splice sites in each CDS at each trinucleotide context
  message("[Step 3/3] Cataloguing the impacts of all possible coding changes in all collected transcripts...")
  
  nt = c("A","C","G","T")
  trinuc_list = paste(rep(nt,each=16,times=1), rep(nt,each=4,times=4), rep(nt,each=1,times=16), sep="")
  trinuc_ind = structure(1:64, names=trinuc_list)
  
  trinuc_subs = NULL; for (j in 1:length(trinuc_list)) { trinuc_subs = c(trinuc_subs, paste(trinuc_list[j], paste(substr(trinuc_list[j],1,1), setdiff(nt,substr(trinuc_list[j],2,2)), substr(trinuc_list[j],3,3), sep=""), sep=">")) }
  trinuc_subsind = structure(1:192, names=trinuc_subs)
  
  # Precalculating a 64x64 matrix with the functional impact of each codon transition (1=Synonymous, 2=Missense, 3=Nonsense)
  impact_matrix = array(NA, dim=c(64,64))
  colnames(impact_matrix) = rownames(impact_matrix) = trinuc_list
  for (j in 1:64) {
    for (h in 1:64) {
      from_aa = seqinr::translate(strsplit(trinuc_list[j],"")[[1]], numcode = numcode)
      to_aa = seqinr::translate(strsplit(trinuc_list[h],"")[[1]], numcode = numcode)
      # Annotating the impact of the mutation
      if (to_aa == from_aa){ 
        impact_matrix[j,h] = 1
      } else if (to_aa == "*"){
        impact_matrix[j,h] = 3
      } else if ((to_aa != "*") & (from_aa != "*") & (to_aa != from_aa)){
        impact_matrix[j,h] = 2
      } else if (from_aa=="*") {
        impact_matrix[j,h] = NA
      }
    }
  }
  
  
  
  build_L_matrix = function(entry) {
    L = array(0, dim=c(192,4))
    cdsseq = as.character(as.vector(entry$seq_cds))
    cdsseq1up = as.character(as.vector(entry$seq_cds1up))
    cdsseq1down = as.character(as.vector(entry$seq_cds1down))
    
    # 1. Exonic mutations
    ind = rep(1:length(cdsseq), each=3)
    old_trinuc = paste(cdsseq1up[ind], cdsseq[ind], cdsseq1down[ind], sep="")
    new_base = c(sapply(cdsseq, function(x) nt[nt!=x]))
    new_trinuc = paste(cdsseq1up[ind], new_base, cdsseq1down[ind], sep="")
    codon_start = rep(seq(1,length(cdsseq),by=3),each=9)
    old_codon = paste(cdsseq[codon_start], cdsseq[codon_start+1], cdsseq[codon_start+2], sep="")
    pos_in_codon = rep(rep(1:3, each=3), length.out=length(old_codon))
    aux = strsplit(old_codon,"")
    new_codon = sapply(1:length(old_codon), function(x) { new_codonx = aux[[x]]; new_codonx[pos_in_codon[x]] = new_base[x]; return(new_codonx) } )
    new_codon = paste(new_codon[1,], new_codon[2,], new_codon[3,], sep="")
    
    imp = impact_matrix[(trinuc_ind[new_codon]-1)*64 + trinuc_ind[old_codon]]
    matrind = trinuc_subsind[paste(old_trinuc, new_trinuc, sep=">")]
    
    # Synonymous
    matrix_ind = table(matrind[which(imp==1)])
    L[as.numeric(names(matrix_ind)), 1] = matrix_ind
    
    # Missense
    matrix_ind = table(matrind[which(imp==2)])
    L[as.numeric(names(matrix_ind)), 2] = matrix_ind
    
    # Nonsense
    matrix_ind = table(matrind[which(imp==3)])
    L[as.numeric(names(matrix_ind)), 3] = matrix_ind
    
    # 2. Splice site mutations
    if (length(entry$intervals_splice)>0) {
      splseq = as.character(as.vector(entry$seq_splice))
      splseq1up = as.character(as.vector(entry$seq_splice1up))
      splseq1down = as.character(as.vector(entry$seq_splice1down))
      old_trinuc = rep(paste(splseq1up, splseq, splseq1down, sep=""), each=3)
      new_trinuc = paste(rep(splseq1up, each=3), c(sapply(splseq, function(x) nt[nt!=x])), rep(splseq1down,each=3), sep="")
      matrind = trinuc_subsind[paste(old_trinuc, new_trinuc, sep=">")]
      matrix_ind = table(matrind)
      L[as.numeric(names(matrix_ind)), 4] = matrix_ind
    }
    
    return(L)
  }
  genes_with_results = gene_list[gene_list %in% names(ref_env)]
  RefCDS = array(list(NULL),length(genes_with_results))
  L_mats = pbapply::pblapply(genes_with_results, function(x) build_L_matrix(ref_env[[x]]), cl = cores)
  for (i in 1:length(genes_with_results)) {
    gene = genes_with_results[i]
    entry = ref_env[[gene]]
    
    RefCDS[[i]] = c(entry, list(L=L_mats[[i]]))
  }
  
  ## Saving the reference GenomicRanges object
  aux = unlist(sapply(1:length(RefCDS), function(x) t(cbind(x,rbind(RefCDS[[x]]$intervals_cds,cbind(RefCDS[[x]]$intervals_splice,RefCDS[[x]]$intervals_splice))))))
  df_genes = as.data.frame(t(array(aux,dim=c(3,length(aux)/3))))
  colnames(df_genes) = c("ind","start","end")
  df_genes$chr = unlist(sapply(1:length(RefCDS), function(x) rep(RefCDS[[x]]$chr,nrow(RefCDS[[x]]$intervals_cds)+length(RefCDS[[x]]$intervals_splice))))
  df_genes$gene = sapply(RefCDS, function(x) x$gene_name)[df_genes$ind]
  
  gr_genes = GenomicRanges::GRanges(df_genes$chr, IRanges::IRanges(df_genes$start, df_genes$end))
  GenomicRanges::mcols(gr_genes)$names = df_genes$gene
  
  return(list(RefCDS, gr_genes))
  
} # EOF
