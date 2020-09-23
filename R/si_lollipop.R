#' Plot selection intensities
#' 
#' Lollipop plots can be used to compare highly selected variants
#' within or among groups of samples.
#' 
#' To keep the lollipops readable, no more than \code{max_sites} variants are shown for
#' any group. At the same time, all SIs that are within the limits of the output plot are
#' depicted (which means that specification of \code{max_sites} creates an implicit SI
#' floor). If the labels don't look good, you may have to reduce max_sites or break 
#' up the plot into pieces. (If you find a better way to get lots of variant labels to
#' display nicely, we would love to hear about it!)
#' 
#' 
#' @import ggplot2
#' @import ggrepel
#' @param si_list Named list of tables, usually created by subsetting from selection
#'   tables. Each table should have "variant_name", "variant_type", "gene", and
#'   "selection_intensity" columns.
#' @param title plot title to pass to ggplot 
#' @param max_sites maximum number of variant sites to include per lollipop;
#'                  if you try to include too many, you may have a challenge
#'                  getting it to look good
#' @return ggplot object with lollipops
#' @examples 
#' \dontrun{
#' # Compare cancer subtypes from different analyses
#' si_lollipops(si_list = list(luad = luad_cesa$selection, lusc = lusc_cesa$selection))
#' 
#' # Compare results for two genes
#' tp53 = cesa$selection[gene == 'TP53']
#' egfr = cesa$selection[gene == 'EGFR']
#' si_lollipops(si_list = list(EGFR = egfr, TP53 = tp53))
#' }
#' @export
si_lollipops = function(si_list, title = "My SIs", max_sites = 50) {
  if (! is(si_list, "list") || ! all(sapply(si_list, is.data.table))) {
    stop("si_list should be a named list of data.tables.")
  }
  num_groups = length(si_list)
  group_names = names(si_list)
  
  if(is.null(group_names)) {
    stop("si_list should be a named list (for labeling lollipops)")
  }
  
  if(! is(title, "character") | length(title) > 1) {
    stop("title should be 1-length character")
  }
  
  if(! is(max_sites, "numeric")) {
    stop("max_sites should be numeric")
  }
  
  # for each group, get the highest N SIs (set by max_sites), then pick the
  # greatest of all of these as the SI floor so that no more than max_sites SIs print per lollipop
  min_si = 0
  for (i in 1:num_groups) {
    curr = si_list[[i]]
    if (curr[, .N] > max_sites) {
      lowest_passing = curr[order(selection_intensity, decreasing = T), selection_intensity[max_sites]]
      if (lowest_passing > min_si) {
        min_si = lowest_passing
      }
    }
  }
  # min_si = max(sapply(si_list, function(x) x[order(selection_intensity, decreasing = T), 
  #                                            min(selection_intensity[1:max_sites], na.rm = T)]))


  placeholder = data.table(group = numeric(), selection_intensity = numeric(), variant_name_print = character(),
                           display_color = factor())
  xlim_low = .8
  xlim_high = num_groups + 1
  lolli = ggplot(placeholder, mapping = aes(group, selection_intensity, label = variant_name_print, color = display_color)) + 
          scale_y_log10() + ylab(expression('Selection intensity'~scriptstyle(~~(log[10])))) +
          theme(
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 16),
            axis.text.x = element_text(size = 16),
            legend.position = "none",
            plot.title = element_text(size = 20, hjust = .5)
          ) + ggtitle(label = title) +  scale_colour_discrete() + 
          scale_x_continuous(breaks = 1:num_groups, labels = group_names, limits = c(xlim_low, xlim_high))
  for (i in 1:num_groups) {
    top_hits = si_list[[i]]
    if (! all(c("selection_intensity", "variant_name", "gene", "variant_type") %in% names(top_hits))) {
      stop("Missing required columns in input table ", i, " (see help).")
    }
    top_hits = top_hits[order(selection_intensity, decreasing = T)][selection_intensity > min_si]
    top_hits[, log_si := log10(selection_intensity)]
    top_hits[, variant_name_print := gsub('_', ' ', variant_name)]
    
    # leave out gene names if all AACs have same gene
    # otherwise, keep gene names and add them to SNVs as well
    if (top_hits[variant_type == "aac", .N] > 0 & top_hits[variant_type == "aac", length(unique(gene)) == 1]) {
      regex = paste0('^', top_hits[variant_type == "aac", gene], ' ')
      old_names = top_hits[variant_type == "aac", variant_name_print]
      new_names = mapply(function(r, n) gsub(r, '', n), regex, old_names)
      top_hits[variant_type == "aac", variant_name_print := new_names]
    } else {
      top_hits[variant_type == "aac", variant_name_print := paste0(' ', variant_name_print)] # to line up with SNV
      # usually there will be an intergenic annotation column
      if ("intergenic" %in% names(top_hits)) {
        top_hits[variant_type == "snv" & intergenic == T, variant_name_print := paste0('(intergenic) ', variant_name_print)]
        top_hits[variant_type == "snv" & intergenic == F, variant_name_print := paste0('(', gene, ') ', variant_name_print)]   
      }
    }
    top_hits[, group := i]
    
    # since intergenic variants have no gene, specify display color
    top_hits[, display_color := factor(gene)]
    
    # shows genes in different colors if more than one
    multiple_genes = length(unique(top_hits$gene)) > 1
    
    # optimize nudge for given number of sites (need more space with more sites)

    if (top_hits[, .N] > 10) {
      # use the plot's "center of mass"
      center = trunc(top_hits[, sum((1:.N) * log10(selection_intensity)) / sum(log10(selection_intensity))])
      top_hits[1:center, nudge := .1 + .55 * (1:.N)/.N]
      top_hits[is.na(nudge), nudge := .1 + .55 * (.N:1)/.N]
    } else {
      top_hits$nudge = .2
    }

    
    #top_hits[, nudge := nudge_x + .2*((-1)^(1:.N))]
    top_hits[, label_size := 3 - .1 * max(nchar(variant_name_print) - 14, 0), by = "variant_name"]
    lolli = lolli + 
            geom_vline(xintercept = i, color = "grey") +
            geom_label_repel(data = top_hits, nudge_x = top_hits$nudge, box.padding = 2, direction = 'y', 
                             segment.size = .3, hjust = 0, color = "black", size = top_hits$label_size) +
            geom_point(data = top_hits, alpha = .7)
    
  }
  return(lolli)
}


