#' cancereffectsizeR's RefCDS builder
#'
#' Based on the buildref function in Inigo Martincorena's package dNdScv, this function
#' takes in gene/transcript/CDS definitions and creates a dNdScv-style RefCDS object and
#' an associated GenomicRanges object also required to run dNdScv.
#'
#' Required columns are seqnames, start, end, strand, gene_name, gene_id, protein_id, and
#' type.  Only rows that have type == "CDS" will be used. Strand should be
#' "+" or "-".
#' 
#' By default, only one the longest complete transcript is used from each gene in the
#' input. If you set use_all_transcripts = TRUE, then all complete transcripts will be
#' used, resulting in multiple RefCDS entries for some genes. If you do this, you may
#' want to first eliminate low-confidence or superfluous transcripts from the input data.
#' 
#' @return A two-item list: RefCDS (which is itself a big list, with each entry containing
#'   information on one coding sequence (CDS)), and a GRanges object that defines the
#'   genomic intervals covered by each CDS.
#' @param gtf Path of a Gencode-style GTF file, or an equivalently formatted data
#'   table. See details for required columns (features). It's possible to build such a 
#'   table using data pulled from biomaRt, but it's easier to use a GTF.
#' @param genome genome assembly name (e.g., "hg19"); an associated BSgenome object must be
#'   available to load.
#' @param use_all_transcripts T/F (default TRUE): Whether to use all complete transcripts or just the longest
#'   one for each gene (defaults to the latter).
#' @param cds_ranges_lack_stop_codons The CDS records in Gencode GTFs don't include the
#'   stop codons in their genomic intervals. If your input does include the stop 
#'   codons within CDS records, set to FALSE.
#' @param cores how many cores to use for parallel computations
#' @param additional_essential_splice_pos Usually not needed. A list of
#'   additional essential splice site positions to combine with those calculated
#'   automatically by this function. Each element of the list should have a name
#'   matching a protein_id in the input and consist of a numeric vector of
#'   additional positions. This option exists so that mutations at chr17:7579312
#'   on TP53 are treated as splice site mutations in cancereffectsizeR's default
#'   hg19 reference data set. (Variants at this coding position, which are
#'   always synonymous, have validated effects on splicing, even though the
#'   position misses automatic "essential splice" annotation by 1 base.)
#' @param numcode (don't use) NCBI genetic code number; currently only code 1, the
#'   standard genetic code, is supported
#' @export

build_RefCDS = function(gtf, genome, use_all_transcripts = TRUE, cds_ranges_lack_stop_codons = TRUE, cores = 1, 
                        additional_essential_splice_pos = NULL, numcode = 1) {
  message("[Step 1/5] Loading data and identifying complete transcripts...")
  
  # Load genome
  bsg = BSgenome::getBSgenome(genome)
  suppressWarnings({GenomeInfoDb::seqlevelsStyle(bsg) = "NCBI"}) # avoid unswitchable seqlevels complaints
  
  # Support for alternate genetic codes not yet implemented
  if(numcode != 1 || length(numcode) != 1) {
    stop("Currently only the standard genetic code is supported.", call. = F)
  }
  
  # validated additional_essential_splice_pos list
  if (! is.null(additional_essential_splice_pos)) {
    if (! is(additional_essential_splice_pos, "list")) {
      stop("additional_essential_splice_pos should be type list")
    }
    tmp = names(additional_essential_splice_pos)
    if(is.null(tmp) || ! length(tmp) == length(unique(tmp))) {
      stop("additional_essential_splice_pos should be a list with uniquely named numeric elements")
    }
    if(! all(sapply(additional_essential_splice_pos, is.numeric))) {
      stop("All elements of additional_essential_splice_pos list should be type numeric.")
    }
  }
  required_cols = c("seqnames", "start", "end", "strand", "gene_id", "gene_name", "protein_id", "type")
  
  
  if (! is.logical(use_all_transcripts) || length(use_all_transcripts) != 1) {
    stop("use_all_transcripts must be TRUE/FALSE.")
  }
  
  if (is(gtf, "character")) {
    if (length(gtf) != 1) {
      stop("gtf should be a file path or data.table", call. = F)
    }
    if (! file.exists(gtf) ) {
      stop("Input GTF file not found.", call. = F)
    }
    gtf = as.data.table(rtracklayer::import(gtf, format = "gtf"))
  } else if (is(gtf, "GRanges")) {
    gtf = as.data.table(gtf)
  } else if (! is(gtf, "data.table")) {
    stop("Input GTF must be text file (GTF format), GRanges, or data.table.")
  }
  if (! all(required_cols %in% names(gtf))) {
    stop("Missing required columns.", call. = F)
  }
  gtf = gtf[, ..required_cols]
  reftable = gtf[type == "CDS", -"type"]
  
  if(reftable[, .N] == 0) {
    stop("No entries of type \"CDS\" in GTF input; these entries are required.")
  }
  
  char_cols = c("seqnames", "gene_id", "gene_name", "strand", "protein_id")
  num_cols = c("start", "end")
  reftable[, (char_cols) := lapply(.SD, as.character), .SDcols = char_cols]
  reftable[, (num_cols) := lapply(.SD, as.numeric), .SDcols = num_cols]
  
  # if strand given as 1/-1, convert to +/-
  reftable[strand == "1", strand := '+']
  reftable[strand == "-1", strand := '-']
  
  # Remove records with chromosomes that don't match reference
  prev_n = reftable[, .N]
  valid_seqnames = as.character(GenomicRanges::seqnames(bsg))
  bad_seq_ind = reftable[ ! seqnames %in% valid_seqnames, which = T]
  if (length(bad_seq_ind) > 0) {
    good_seq = reftable[! bad_seq_ind]
    bad_seq = reftable[bad_seq_ind]
    # Try to help out by stripping chr prefixes and dealing with chrM vs MT on hg19/hg38
    if (genome %in% c("hg19", "hg38")) {
      bad_seq[seqnames == "chrM", seqnames := "MT"]
    }
    bad_seq[, seqnames := gsub("^chr", "", seqnames)]
    still_bad_ind = bad_seq[! seqnames %in% valid_seqnames, which = T]
    now_good_seq = bad_seq[! still_bad_ind]
    reftable = rbind(good_seq, now_good_seq)
    num_still_bad = length(still_bad_ind)
    if (num_still_bad > 0) {
      message(sprintf("Excluded %i records for having seqnames not present in reference:", num_still_bad))
      bad_chr = bad_seq[still_bad_ind, unique(seqnames)]
      message("\t", paste(bad_chr, collapse = ", "))
    }
  }
  
  # Gencode GTFs give out-of-chromosome coordinates on Y-chromosome pseudoautosomal genes
  if (genome %in% c("hg19", "hg38")) {
    par_gene_ind = reftable[grepl('_PAR_Y$', gene_id), which = T]
    num_par = length(par_gene_ind)
    if (num_par > 0) {
      reftable = reftable[!par_gene_ind]
      message(num_par, " pseudoautosomal genes associated with chrY were dropped in favor of their chrX counterparts.")
    }
  }
  
  if (anyNA(reftable)) {
    stop("Transcript data has NA entries in required fields of CDS entries.", call. = F)
  }
  
  # order by protein ID, then exon coding order so that all exons are in transcription order regardless of strand
  reftable = reftable[order(start)]
  reftable[strand == '+', exon_order :=  seq_len(.N), by = "protein_id"]
  reftable[strand == "-", exon_order := rev(seq_len(.N)), by = "protein_id"]
  reftable = reftable[order(protein_id, exon_order)]
  reftable[, width := end - start + 1]
  reftable[, cds_start := cumsum(width) - width + 1, by = "protein_id"]
  reftable[, cds_end := cds_start + width - 1]
  reftable[, cds_length := sum(width), by = "protein_id"]
  
  
  # Gencode GTFs (and possibly others) don't include the stop codon in CDS sequences. 
  # We need these, so we will add them in and edit ranges when cds_ranges_lack_stop_codons = T.
  # But first, we'll check 100 random final exons from decently-sized transcripts and see if the
  # user's designation seems correct.
  to_check = reftable[cds_end == cds_length & cds_length > 1000][sample(1:.N, min(100, .N))]
  
  # Edge case: Sometimes people will run this on a single gene (or set of genes) shorter than 1000
  if (to_check[, .N] == 0) {
    to_check = reftable[cds_end == cds_length][sample(1:.N, min(100, .N))]
  }
  
  to_check[strand == "+", c("codon_start", "codon_end") := list(end - 2, end)]
  to_check[strand == "-", c("codon_start", "codon_end") := list(start, start + 2)]
  final_codon_gr = GenomicRanges::makeGRangesFromDataFrame(to_check[, .(chr = seqnames, start = codon_start, end = codon_end, strand)],
                                                           seqinfo = seqinfo(bsg))
  seqs_to_check = BSgenome::getSeq(bsg, final_codon_gr)
  last_codons = Biostrings::translate(seqs_to_check, no.init.codon = T, if.fuzzy.codon = "X") # avoid errors on Ns
  frac_stop = mean(last_codons == "*")
  if (frac_stop > .5 & cds_ranges_lack_stop_codons == TRUE) {
    stop(paste0("You ran with cds_ranges_lack_stop_codons = T, but a quick check of random transcripts found\n",
                   sprintf("that %.0f%% of them had terminating stop codons.", frac_stop * 100)))
  } else if (frac_stop < .5 & cds_ranges_lack_stop_codons == FALSE) {
    msg = paste0("You ran with cds_ranges_lack_stop_codons = F, but a quick check of random transcripts found\n",
                sprintf("that only %.0f%% of them had terminating stop codons.", frac_stop * 100))
    if (frac_stop < .1) {
      stop(msg)
    }
    warning(msg, immediate. = T)
  }
  
  # Expand ranges to include stop codons if needed (which also means adding 3 to length)
  if (cds_ranges_lack_stop_codons) {
    reftable[, cds_length := cds_length + 3]
    reftable[cds_length - cds_end == 3, cds_end := cds_length]
    reftable[cds_end == cds_length & strand == "-", start := start - 3]
    reftable[cds_end == cds_length & strand == "+", end := end + 3]
  }
  
  # change a few column names for clarity
  setnames(reftable, c("seqnames", "start", "end"), c("chr", "genomic_start", "genomic_end"))
  
  ## dNdScv's buildref function gives a warning if some gene names map to more than one
  ## gene ID, and offers a "longname" option When using a Gencode GTF, PAR genes get
  ## multiple gene IDs (one for the chrX version, one for the chrY), and this script
  ## already handles those by taking the chrX version only. Beyond that, recent Gencode
  ## GTFs have very few gene names with multiple gene IDs (like, around 10), so we're not
  ## going to worry about them.
  
  # longname = paste(reftable$gene.id, reftable$gene_name, sep=":") # Gene name combining the Gene stable ID and the Associated gene name
  # if (useids==T) {
  #   reftable$gene_name = longname # Replacing gene names by a concatenation of gene ID and gene name
  # }
  # if (length(unique(reftable$gene_name))<length(unique(longname))) {
  #   warning(sprintf("%0.0f unique gene IDs (column 1) found. %0.0f unique gene names (column 2) found. Consider combining gene names and gene IDs or replacing gene names by gene IDs to avoid losing genes (see useids argument in ? buildref)",length(unique(reftable$gene.id)),length(unique(reftable$gene_name))))
  # }
  # 
  
  # Confusing, but for transcript-level compatibility with dNdScv, which finds CDS entries
  # by the "gene_name" attribute, we need to call the transcript IDs gene_name in the
  # output, and add an extra attribute with the real gene name
  if (use_all_transcripts) {
    reftable[, entry_name := protein_id]
  } else {
    reftable[, entry_name := gene_name]
  }
  
  # Removing CDS of length not multiple of 3, and sort CDS from longest to shortest for each feature
  # To-do: Investigate the CDS sequences that don't pass this
  reftable = reftable[cds_length %% 3 == 0]
  reftable = reftable[order(gene_name, -cds_length, exon_order)]
  all_entries = unique(reftable$entry_name)
  
  message("[Step 2/5] Processing exons...")
  grs = GenomicRanges::makeGRangesFromDataFrame(reftable, keep.extra.columns = T, seqinfo = GenomicRanges::seqinfo(bsg), 
                                                seqnames.field = "chr", start.field = "genomic_start", end.field = "genomic_end", 
                                                strand.field = "strand")
  GenomicRanges::start(grs) = GenomicRanges::start(grs) - 1
  GenomicRanges::end(grs) = GenomicRanges::end(grs) + 1
  padded_seqs = BSgenome::getSeq(bsg, grs)
  padded_cds_seqs = S4Vectors::split(padded_seqs, f = reftable$protein_id)
  
  remaining_entries = all_entries
  curr_table = reftable
  
  final_cdsseq = DNAStringSet()
  final_cdsseq1up = DNAStringSet()
  final_cdsseq1down = DNAStringSet()
  
  while(length(remaining_entries) > 0) {
    setkey(curr_table, "entry_name")
    curr_cds = curr_table[remaining_entries, , mult = "first"]
    setkey(curr_cds, "entry_name")

    # CDS sequences have 1-base padding on start and end of each interval range
    curr_cds_seqs = padded_cds_seqs[curr_cds$protein_id]
    
    
    all_seqs_as_list = S4Vectors::lapply(curr_cds_seqs, function(x) list(unlist(Biostrings::subseq(x, 2, -2)),
                                                                         unlist(Biostrings::subseq(x, 1, -3)),
                                                                         unlist(Biostrings::subseq(x, 3, -1))))
    
    cdsseq = DNAStringSet(lapply(all_seqs_as_list, function(x) x[[1]]))
    
    # Note that we're not checking the upstream/downstream bases, which could have N's
    # Since the original buildref didn't account for this possibility, we won't either, yet
    cdsseq = cdsseq[Biostrings::vcountPattern("N", cdsseq) == 0]
    pseq = Biostrings::translate(cdsseq)
    
    # we don't want any stop codons before the last elements of the AAStrings
    good_cds = names(pseq[Biostrings::vcountPattern("*", Biostrings::subseq(pseq, 1, -2)) == 0])
    
    cdsseq = cdsseq[good_cds]
    cdsseq1up = DNAStringSet(lapply(all_seqs_as_list[good_cds], function(x) x[[2]]))
    cdsseq1down = DNAStringSet(lapply(all_seqs_as_list[good_cds], function(x) x[[3]]))

    final_cdsseq = c(final_cdsseq, cdsseq)
    final_cdsseq1up = c(final_cdsseq1up, cdsseq1up)
    final_cdsseq1down =  c(final_cdsseq1down, cdsseq1down)
    
    # remove entries used on last attempt from the table
    curr_table = curr_table[! curr_cds$protein_id, on = "protein_id"]
    
    # When doing one CDS per gene, drop remaining entries after getting a good one in a gene
    if (! use_all_transcripts) {
      remaining_entries = setdiff(remaining_entries, curr_cds[names(cdsseq), entry_name, on = 'protein_id'])
    }
    # Some genes will probably not be found (i.e., no CDS was good for them)
    # Take intersection with remaining genes/proteins in the table to get those that could still have a good CDS
    remaining_entries = intersect(remaining_entries, curr_table$entry_name)
    
    # For efficiency, drop records for any genes not yet to be completed
    curr_table = curr_table[remaining_entries, on = 'entry_name']
  }
  
  # Convert to environments, which are hashed, for faster retrieval
  final_cdsseq = list2env(as.list(final_cdsseq))
  final_cdsseq1up = list2env(as.list(final_cdsseq1up))
  final_cdsseq1down =  list2env(as.list(final_cdsseq1down))
  
  message("[Step 3/5] Handling splice sites...")
  # Keep just the reftable entries of the accepted CDS
  reftable = reftable[names(final_cdsseq), on = "protein_id"]
  setkey(reftable, "protein_id")
  
  # Get splice site info for the final CDS set
  # Get essential splice sites: 1,2 bp upstream of exon start (5') and 1,2,5 bp downstream of exon end (3')
  # These are described as the most important positions for correct splicing in 10.1534/genetics.105.044677,
  # but no reasoning or citations are given.
  # Note that single-exon transcripts have no splice sites; they'll be processed incorrectly and then fixed below
  spl_neg = data.table()
  if(reftable[strand == '-', .N] > 0) {
    spl_neg = reftable[strand == '-', .(chr = chr[1], strand = strand[1], us_ends = list(genomic_end[2:.N]), ds_ends = list(genomic_start[1:(.N - 1)])), by = "protein_id"]
    spl_neg[, splpos := list(list(sort(unique(c(unlist(us_ends) + 1, unlist(us_ends) + 2, unlist(ds_ends) - 1, unlist(ds_ends) - 2, unlist(ds_ends) - 5))))),
            by = "protein_id"]
  }
  spl_positive = data.table()
  if(reftable[strand == '+', .N] > 0) {
    spl_positive = reftable[strand == '+', .(chr = chr[1], strand = strand[1], us_ends = list(genomic_start[2:.N]), ds_ends = list(genomic_end[1:(.N - 1)])), by = "protein_id"]
    spl_positive[, splpos := list(list(sort(unique(c(unlist(us_ends) - 1, unlist(us_ends) - 2, unlist(ds_ends) + 1, unlist(ds_ends) + 2, unlist(ds_ends) + 5))))),
                 by = "protein_id"]
  }
  spl = rbind(spl_positive, spl_neg)
  setkey(spl, "protein_id")
  
  # Insert user-supplied essential splice positions  in addition to those inferred from distance to splice site
  if (! is.null(additional_essential_splice_pos)) {
    extra_pid = names(additional_essential_splice_pos)
    for (i in 1:length(extra_pid)){
      curr_pid = extra_pid[[i]]
      prev_splpos = spl[curr_pid, unlist(splpos)]
      if (is.null(prev_splpos)) {
        warning("Additional essential splice sites supplied for ", curr_pid, " but this protein ID does not appear in output.")
        next
      }
      new_splpos = additional_essential_splice_pos[[curr_pid]]
      updated_splpos = sort(unique(c(prev_splpos, new_splpos)))
      spl[curr_pid, splpos := list(updated_splpos)]
    }
  }
  
  
  # remove single-exon records, get splice site squences, and save splice site positions
  single_exon_pids = reftable[, .N, by = "protein_id"][N == 1, protein_id]
  spl = spl[!single_exon_pids]
  
  for_spl_gr = spl[, .(chr, strand, pos = unlist(splpos)), by = "protein_id"]
  for_spl_gr[, c("start", "end") := .(pos - 1, pos + 1)]
  spl_gr = GenomicRanges::makeGRangesFromDataFrame(for_spl_gr, strand.field = "strand", seqinfo = GenomicRanges::seqinfo(bsg))
  
  padded_spl_seqs = BSgenome::getSeq(bsg, spl_gr)
  padded_spl_seqs = S4Vectors::split(padded_spl_seqs, f = for_spl_gr$protein_id)
  rm(for_spl_gr)
  
  all_seqs_as_list = S4Vectors::lapply(padded_spl_seqs, function(x) list(unlist(Biostrings::subseq(x, 2, -2)),
                                                                       unlist(Biostrings::subseq(x, 1, -3)),
                                                                       unlist(Biostrings::subseq(x, 3, -1))))
  
  splseq = list2env(lapply(all_seqs_as_list, function(x) x[[1]]))
  splseq1up = list2env(lapply(all_seqs_as_list, function(x) x[[2]]))
  splseq1down = list2env(lapply(all_seqs_as_list, function(x) x[[3]]))
  
  spl_names = spl$protein_id
  spl = spl$splpos
  names(spl) = spl_names
  
  message("[Step 4/5] Initializing RefCDS object...")
  
  # Gather info for each gene (or transcript) and convert strand from +/- to numeric 1/-1
  # Also hold grab the real gene name when "gene_name" actually refers to transcript ID
  gene_info = reftable[names(final_cdsseq), .(gene_name, gene_id, protein_id, CDS_length = cds_length, chr, 
                                                char_strand = as.character(strand), entry_name), mult = "first"]
  
  # RefCDS entries must match format expected by dNdScv (see message at the end)
  if (use_all_transcripts) {
    gene_info[, c("real_gene_name", "gene_name") := list(gene_name, entry_name)]
  }

  gene_info[char_strand == "-", strand := -1]
  gene_info[char_strand == "+", strand := 1]
  gene_info[, char_strand := NULL]
  setorder(gene_info, entry_name) # reorder rows by gene name
  
  # to meet RefCDS spec, re-order table so that CDS intervals are in genomic rather than coding order
  reftable = reftable[order(protein_id, genomic_start)]
  cds_intervals = split(reftable[, .(genomic_start, genomic_end)], f = reftable$protein_id)
  cds_intervals = lapply(cds_intervals, function(x) unname(as.matrix(x)))

  refcds = list()
  for (i in 1:nrow(gene_info)) {
    current = as.list(gene_info[i])
    entry = current$entry_name
    current[["entry_name"]] = NULL
    pid = current$protein_id
    current[["intervals_cds"]] = cds_intervals[[pid]]
    
    intervals_splice = spl[[pid]]
    has_splice = TRUE
    if (is.null(intervals_splice)) {
      current[["intervals_splice"]] = numeric()
      no_splice = FALSE
    } else {
      current[["intervals_splice"]] = intervals_splice
    }
    
    current[["seq_cds"]] = final_cdsseq[[pid]]
    current[["seq_cds1up"]] = final_cdsseq1up[[pid]]
    current[["seq_cds1down"]] = final_cdsseq1down[[pid]]
    
    if (has_splice) {
      current[["seq_splice"]] = splseq[[pid]]
      current[["seq_splice1up"]] = splseq1up[[pid]]
      current[["seq_splice1down"]] = splseq1down[[pid]]
    }
    
    refcds[[entry]] = current
  }

  ## L matrices: number of synonymous, missense, nonsense and splice sites in each CDS at each trinucleotide context
  message("[Step 5/5] Cataloguing the impacts of all possible coding changes in all collected transcripts...")
  
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
  
  other_nt = list(A = c("C", "G", "T"), C = c("A", "G", "T"), G = c("A", "C", "T"), 
             T = c("A", "C", "G"))
  build_L_matrix = function(entry) {
    L = array(0, dim=c(192,4))
    cdsseq = as.character(as.vector(entry$seq_cds))
    cdsseq1up = as.character(as.vector(entry$seq_cds1up))
    cdsseq1down = as.character(as.vector(entry$seq_cds1down))
    
    # 1. Exonic mutations
    ind = rep(1:length(cdsseq), each=3)
    old_trinuc = paste(cdsseq1up[ind], cdsseq[ind], cdsseq1down[ind], sep="")
    new_base = c(sapply(cdsseq, function(x) other_nt[[x]]))
    new_trinuc = paste(cdsseq1up[ind], new_base, cdsseq1down[ind], sep="")
    codon_start = rep(seq(1,length(cdsseq),by=3),each=9)
    old_codon = paste(cdsseq[codon_start], cdsseq[codon_start+1], cdsseq[codon_start+2], sep="")
    pos_in_codon = rep(rep(1:3, each=3), length.out=length(old_codon))
    new_codon = old_codon
    substr(new_codon, pos_in_codon, pos_in_codon) = new_base

    
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
  entries_with_results = names(refcds)
  L_mats = pbapply::pblapply(entries_with_results, function(x) build_L_matrix(refcds[[x]]), cl = cores)
  for (i in 1:length(entries_with_results)) {
    entry = entries_with_results[i]
    refcds[[entry]][["L"]] = L_mats[[i]]
  }
  refcds = as.array(refcds) # historically motivated conversion
  
  ## Saving the reference GenomicRanges object
  aux = unlist(sapply(1:length(refcds), function(x) t(cbind(x,rbind(refcds[[x]]$intervals_cds,cbind(refcds[[x]]$intervals_splice,refcds[[x]]$intervals_splice))))))
  df_genes = as.data.frame(t(array(aux,dim=c(3,length(aux)/3))))
  colnames(df_genes) = c("ind","start","end")
  df_genes$chr = unlist(sapply(1:length(refcds), function(x) rep(refcds[[x]]$chr,nrow(refcds[[x]]$intervals_cds)+length(refcds[[x]]$intervals_splice))))
  df_genes$entry = names(refcds)[df_genes$ind]
  
  # note that gr_cds is what should be supplied to dNdScv as "gr_genes"
  gr_cds = GenomicRanges::GRanges(df_genes$chr, IRanges::IRanges(df_genes$start, df_genes$end))
  GenomicRanges::mcols(gr_cds)$names = df_genes$entry
  
  if (use_all_transcripts) {
    # since df_genes$gene is actually a protein ID, match with the "gene_name" in 
    # df_genes (also protein ID), and get corresponding real_gene_name
    actual_gene_names = gene_info[df_genes, real_gene_name, on = c(entry_name = "entry")]
    GenomicRanges::mcols(gr_cds)$gene = actual_gene_names
    pretty_message(paste0("Note: For compatibility with dNdScv, the \"gene_name\" attribute of each RefCDS entry is ",
                   "actually the protein ID (sorry). A \"real_gene_name\" attribute has also been added ",
                   "to each entry, and cancereffectsizeR will automatically keep things straight."))
  }
  return(list(refcds, gr_cds))
}
