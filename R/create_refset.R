#' Create a custom refset
#'
#' Use this function to create and save a directory of custom reference data that can be
#' used with cancereffectsizeR instead of supplied refsets like \code{ces.refset.hg19}. All
#' arguments are required except default_exome/exome_interval_padding, which are recommended.
#' 
#' To run this function, you'll need to have output from \code{build_RefCDS()}.
#' 
#' @param output_dir Name/path of an existing, writable output directory where all data
#'   will be saved. The name of this directory will serve as the name of the custom refset.
#' @param refcds_dndscv Transcript information in the two-item list (consisting of RefCDS
#'   and gr_genes) that is output by \code{build_RefCDS}. This transcript information will be used with dNdScv.
#' @param refcds_anno Transcript information in the two-item list (consisting of RefCDS
#'   and gr_genes) that is output by \code{build_RefCDS}. This transcript information will be used for
#'   cancereffectsizeR's annotations. If unspecified, the same reference information as supplied for dNdScv will be used.
#' @param species_name Name of the species, primarily for display (e.g., "human").
#' @param genome_build_name Name of the genome build, primarily for display (e.g., "hg19").
#' @param BSgenome_name The name of the BSgenome package to use (e.g., "hg19"); will used by
#'   cancereffectsizeR to load the reference genome via BSgenome::getBSgenome().
#' @param supported_chr Character vector of supported chromosomes. Note that cancereffectsizeR uses
#'   NCBI-style chromosome names, which means no chr prefixes ("X", not "chrX"). Mitochondrial
#'   contigs shouldn't be included since they would require special handling that hasn't been
#'   implemented.
#' @param default_exome A BED file or GRanges object that defines coding regions in the genome as
#'   might be used by an exome capture kit. This file (or GRanges) might be acquired or generated
#'   from exome capture kit documentation, or alternatively, coding regions defined in a GTF file
#'   (or the granges output by build_RefCDS()).
#' @param exome_interval_padding Number of bases to pad start/end of each covered interval, to allow
#'   for some variants to be called just outside of targeted regions, where there still may be
#'   pretty good sequencing coverage.
#' @param transcript_info Additional information about coding (and, optionally, noncoding)
#'   transcripts from a Gencode GTF, supplied as a data.table. See the format provided in
#'   ces.refset.hg38. You'll have to match the format (including column names) pretty closely to get
#'   expected behavior. Noncoding transcripts are represented only by records with transcript_type =
#'   "transcript", and protein-coding transcripts are representing with transcript, CDS, and UTR
#'   records. Note that in Gencode format.
#' @param cores How many cores to use (default 1).
#' @export
create_refset = function(output_dir, refcds_dndscv, refcds_anno = NULL, species_name, genome_build_name, 
                         BSgenome_name, supported_chr = c(1:22, 'X', 'Y'), default_exome = NULL,
                         exome_interval_padding = 0, transcript_info = NULL, cores = 1) {
  if(! is.numeric(cores) || cores - as.integer(cores) != 0) {
    stop('cores should be type integer')
  }
  cores_available = parallel::detectCores()
  if(is.na(cores_available)) cores_available = 1
  if(cores > cores_available) {
    message('Running with cores = ', cores_available, ' since that is how many appear available.')
  }
  
  if(! is.null(transcript_info)) {
    if(! is.data.table(transcript_info)) {
      stop('transcript_info must be a data.table (and a quite specifically formatted one at that).')
    }
  }
  if(! is(refcds_dndscv, 'list') || length(refcds_dndscv) != 2) {
    stop('refcds_dndscv does not look right. Expected 2-length list (containing RefCDS and gr_genes).')
  }
  use_separate_refcds_anno = TRUE
  if(is.null(refcds_anno)) {
    refcds_anno = refcds_dndscv
    use_separate_refcds_anno = FALSE
  } else  if(! is(refcds_anno, 'list') || length(refcds_dndscv) != 2) {
    stop('refcds_anno does not look right. Expected 2-length list (containing RefCDS and gr_genes).')
  }
  if (! is.character(output_dir) || length(output_dir) != 1) {
    stop("output_dir should be 1-length character (path for a new directory)")
  }
  
  if(! dir.exists(output_dir)) {
    stop("The supplied output_dir doesn't exist.")
  }
  
  if(file.access(output_dir, mode = 2) != 0) {
    stop("output_dir doesn't appear to be writable.")
  }
  
  if (! is.character(species_name) || length(species_name) != 1) {
    stop("species_name should be 1-length character")
  }
  if (! is.character(genome_build_name) || length(genome_build_name) != 1) {
    stop("genome_build_name should be 1-length character")
  }
  
  if (! is.character(BSgenome_name) || length(BSgenome_name) != 1) {
    stop("BSgenome_name should be 1-length character")
  }
  
  if (! is.character(supported_chr) || length(supported_chr) == 0) {
    stop("supported_chr should be type character")
  }
  if(is.null(default_exome)) {
    warning("No default exome supplied (okay, but recommended to include if available).")
  } else if(is.character(default_exome)) {
    if(length(default_exome) != 1) {
      stop("default_exome should be 1-length character (file path).")
    }
    if (! file.exists(default_exome)) {
      stop("File specified in default_exome doesn't exist.")
    }
  } else if (! is(default_exome, "GRanges")) {
    stop("default_exome should be a BED file path or GRanges object.")
  }
  
  if (! is.numeric(exome_interval_padding) || length(exome_interval_padding) != 1) {
    stop("exome_interval_padding should be 1-length numeric")
  }
  if (exome_interval_padding > 1000) {
    warning("You chose a high exome_interval_padding!")
  }
  
  # save an environment (dictionary) with some basic genome info
  genome_info = new.env(parent = emptyenv()) # emptyenv or when you call saveRDS, parents will get sucked in
  
  # for displaying to user
  genome_info[['build_name']] = genome_build_name
  genome_info[['species']] = species_name
  
  # name to use when loading genome with BSgenome::getBSgenome()
  genome_info[['BSgenome']] = BSgenome_name
  genome_info[['supported_chr']] = supported_chr
  
  # use seqlevelsStyle from RefCDS, since it will need to match
  chromosome_style = seqlevelsStyle(refcds_anno[[2]])[1]
  genome_info[['chromosome_style']] = chromosome_style
  
  # Compute and save genome-wide trinucleotide contexts
  # (For ces.refset.hg19, this exactly reproduces the existing hg19 counts packaged with deconstructSigs.)
  bsg = BSgenome::getBSgenome(genome_info[['BSgenome']])
  withCallingHandlers(
    { 
      GenomeInfoDb::seqlevelsStyle(bsg) = chromosome_style
    },
    warning = function(w) {
      if (grepl("more than one best sequence renaming map", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      } else if(grepl("cannot switch some of.*to .*style", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }        
  )
  
  missing_chr = setdiff(supported_chr, seqnames(bsg))
  if (length(missing_chr) > 0) {
    stop('Not all supported_chr are present in genome. Check that chromosome naming style is consistent.\n',
         'Missing: ', paste(missing_chr, collapse = ", "), '.')
  }
  
  # deconstructSigs_trinuc_string is internal to cancereffectsizeR (find with getAnywhere if necessary).
  # Here, we need unique trinucleotide contexts (without mutations) in deconstructSigs ordering,
  # which is why we produce this by converting from their column headings.
  # Convert stuff like "A[C>G]A" format to "ACA" (just the reference trinucleotide context).
  context_names = unique(sub("\\[([ACTG]).*\\]", "\\1", deconstructSigs_trinuc_string))
  context_names = sort(context_names)
  reverse_complement_names = unique(as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(context_names))))
  
  message("Loading reference genome and counting trinucleotide contexts...")
  wg = BSgenome::getSeq(bsg)
  genome_tri = Biostrings::trinucleotideFrequency(wg[supported_chr])
  genome_tri = colSums(genome_tri)
  
  genome_counts = genome_tri[context_names] + genome_tri[reverse_complement_names]
  genome_counts = data.frame(x = genome_counts)
  
  
  # Create trinuc context counts for CDS sequences
  # There will be redundant counts for each possible substitution in a given context
  trinuc_mut_names = sub("\\[([ACTG]).*\\]", "\\1", deconstructSigs_trinuc_string)
  
  # deconstructSigs only includes C/T as central nucleotides; we need reverse complement for A/G in center
  reverse_trinuc_mut_names = as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(trinuc_mut_names)))
  
  
  # Get the reference dinucleotides for the COSMIC DBS classes, in order
  dn_names = sub(">.*", "", cosmic_dbs_classes)
  cosmic_dbs_dinuc = unique(dn_names)
  reverse_dn_names = as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(dn_names)))
  
  # We use the annotation RefCDS, not the dNdScv one, because we need context counts for every coding annotation
  message("Counting trinculeotide contexts in ", format(length(refcds_anno[[1]]), big.mark = ','), " RefCDS entries...")
    
  
  process_cds = function(entry) {
    total_counts = integer(96) # number of distinct deconstructSigs sbs
    dn_counts = integer(length(cosmic_dbs_classes))
    # for each transcripts, need to consider each exon and 1 base up/downstream for trinucleotide context
    intervals = entry$intervals_cds # two-column matrix with start/stop coordinates of each cds
    cds_lengths = intervals[,2] - intervals[,1] 
    
    # if transcript is on negative strand, flip exon order
    if (entry$strand == -1) {
      cds_lengths = rev(cds_lengths)
    }
    start = 1
    for (cds_length in cds_lengths) {
      end = start + cds_length
      # xscat and subseq are much more efficient than plain concatenation and subsetting
      seq = Biostrings::xscat(Biostrings::subseq(entry$seq_cds1up, start = start, width = 1), 
                              Biostrings::subseq(entry$seq_cds, start = start, end = end), 
                              Biostrings::subseq(entry$seq_cds1down, start = end, width = 1))
      
      # get trinucleotide counts for cds
      # this function returns full 64-cell table, including counts of 0 (otherwise subsetting below wouldn't work)
      cds_counts = Biostrings::trinucleotideFrequency(seq)
      
      # order the counts in deconstructSigs order and add them to total counts for this transcript
      total_counts = cds_counts[trinuc_mut_names] + cds_counts[reverse_trinuc_mut_names] + total_counts
      
      # Repeat for dinucleotides
      dn_seq = Biostrings::xscat(Biostrings::subseq(entry$seq_cds1up, start = start, width = 1), 
                                 Biostrings::subseq(entry$seq_cds, start = start, end = end), 
                                 Biostrings::subseq(entry$seq_cds1down, start = end, width = 1))
      curr_dn_counts = Biostrings::dinucleotideFrequency(dn_seq)
      
      cosmic_dinuc_counts = curr_dn_counts[cosmic_dbs_dinuc]
      other_dinuc_names = setdiff(names(curr_dn_counts), cosmic_dbs_dinuc)
      other_dinuc_counts = curr_dn_counts[other_dinuc_names]
      other_names_reversed = as.character(reverseComplement(DNAStringSet(other_dinuc_names)))
      cosmic_dinuc_counts[other_names_reversed] = cosmic_dinuc_counts[other_names_reversed] + other_dinuc_counts
      
      dn_counts = dn_counts + cosmic_dinuc_counts[dn_names]
      
      
      # next cds sequence starts with next base in the seq_cds sequence
      start = end + 1
    }
    
    # Need nonzero counts for every context in case the gene's rate is used for a nearby noncoding site
    # with a context that doesn't appear in the CDS.
    if(0 %in% total_counts) {
      total_counts = total_counts + 1
    }
    if(0 %in% dn_counts) {
      dn_counts = dn_counts + 1
    }
    
    # normalize and add trinuc comp to environment (dropping names since they don't match deconstructSigs format)
    comp = unname(total_counts / sum(total_counts))
    dn_comp = unname(dn_counts / sum(dn_counts))
    
    return(list(comp, dn_comp))
  }
  
  trinuc_and_dinuc = pbapply::pblapply(refcds_anno[[1]], process_cds, cl = cores)
  
  
  gene_trinuc_comp = lapply(trinuc_and_dinuc, '[[', 1)
  names(gene_trinuc_comp) = names(refcds_anno[[1]])
  gene_trinuc_comp = list2env(gene_trinuc_comp, parent = emptyenv())
  
  dbs_comp = lapply(trinuc_and_dinuc, '[[', 2)
  names(dbs_comp) = names(refcds_anno[[1]])
  dbs_comp = list2env(dbs_comp, parent = emptyenv())
  
  
  
  # Optional: load and save default exome intervals
  # If you don't set a default exome, the user must always supply coverage intervals.
  using_exome_files = ! is.null(default_exome)
  if (using_exome_files) {
    message("Loading default exome and counting trinucleotide contexts...")
    
    if(is.character(default_exome)) {
      default_exome = rtracklayer::import(default_exome, format = "bed")
    }
    
    genome(default_exome) = genome(bsg)[1]
    suppressWarnings({seqlevelsStyle(default_exome) = seqlevelsStyle(bsg)})
    default_exome = default_exome[as.character(seqnames(default_exome)) %in% supported_chr]
    default_exome = unstrand(default_exome)
    # need to reorder seqlevels (1, 2, ..., X, Y) in cases where the input bed file wasn't sorted
    seqlevels(default_exome) = supported_chr
    default_exome = reduce(sort(default_exome), drop.empty.ranges = T)
    
    # Add interval padding (often 100bp)
    start(default_exome) = start(default_exome) - exome_interval_padding
    end(default_exome) = end(default_exome) + exome_interval_padding
    default_exome = reduce(default_exome, drop.empty.ranges = T)
    
    seqinfo(default_exome) = seqinfo(bsg)
    default_exome = keepSeqlevels(default_exome, supported_chr) # drop extra seqlevels
    
    # get generic captured exome sequence and tabulate trinuc contexts
    exome_seq = getSeq(bsg, default_exome)
    exome_tri_contexts = Biostrings::trinucleotideFrequency(exome_seq)
    exome_tri_contexts = colSums(exome_tri_contexts)
    
    # reorder the counts as desired, then save as a data frame since that's what deconstructSigs wants
    exome_counts = exome_tri_contexts[context_names] + exome_tri_contexts[reverse_complement_names]
    exome_counts = data.frame(x = exome_counts)
  }
  
  gr_genes = refcds_dndscv[[2]]
  if (is.null(gr_genes$gene)) {
    gene_names = unique(gr_genes$names)
  } else {
    gene_names = unique(gr_genes$gene)
  }
  
  if(use_separate_refcds_anno) {
    gr_genes_anno = refcds_anno[[2]]
    if (is.null(gr_genes_anno$gene)) {
      anno_gene_names = sort(unique(gr_genes_anno$names))
    } else {
      anno_gene_names = sort(unique(gr_genes_anno$gene))
    }
    if(! identical(anno_gene_names, sort(gene_names))) {
      stop('The separate RefCDS that were supplied for annotation and dNdScv do not have exactly the same gene names.')
    }
    # If refcds_anno is not the same object as refcds_dndscv, we'll get rid of elements not needed for our annotation.
    refcds_anno[[1]] = lapply(refcds_anno[[1]], '[', c('gene_name', 'gene_id', 'protein_id', 'CDS_length', 'chr', 
                                                       'real_gene_name', 'strand', 'intervals_cds', 'intervals_splice', 'seq_cds'))
  }
  
  message("Saving data files...")
  output_dir = sub('/+$', '', output_dir)
  
  if(using_exome_files) {
    saveRDS(default_exome, paste0(output_dir, "/default_exome_gr.rds"))
    saveRDS(exome_counts, paste0(output_dir, "/tri.counts.exome.rds"))
  }
  
  saveRDS(genome_info , paste0(output_dir, "/genome_build_info.rds"))
  saveRDS(genome_counts, paste0(output_dir, "/tri.counts.genome.rds"))
  saveRDS(refcds_anno[[1]], paste0(output_dir, "/RefCDS.rds"))
  saveRDS(refcds_anno[[2]], paste0(output_dir, "/gr_genes.rds"))
  
  # If we are using separate RefCDS for cancereffectsizeR annotation and dNdScv,
  # that means we still need to save copies of the dNdScv versions.
  if(use_separate_refcds_anno) {
    saveRDS(refcds_dndscv[[1]], paste0(output_dir, "/RefCDS.dndscv.rds"))
    saveRDS(refcds_dndscv[[2]], paste0(output_dir, "/gr_genes.dndscv.rds"))
  }
  
  saveRDS(gene_names, paste0(output_dir, "/gene_names.rds"))
  saveRDS(gene_trinuc_comp, paste0(output_dir, "/gene_trinuc_comp.rds"))
  saveRDS(dbs_comp, paste0(output_dir, '/cds_dbs_exposure.rds'))
  
  if(! is.null(transcript_info)) {
    saveRDS(transcript_info, paste0(output_dir, '/transcript_info.rds'))
  }
}

