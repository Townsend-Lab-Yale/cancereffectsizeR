# File: plot_signature_effects.R
# Author: Derek Song
# Date: August 14, 2023 (last updated 2.26.24 by DS)
# Purpose: Plot and compare the source weight and effect share of mutational
# signatures in one or many groups. To be used with mutational_signature_effects().

# load graphing libraries
library(viridis)
library(ggplot2)

#' @param mutational_effects
#' The output of mutational_signature_effects().
#' Calculate signature source weights and effect shares for each group,
#' then pass all results to the function in either a named list (if there are 
#' multiple groups), or the normal output of mutational_signature_effects
#' (if there is only one group).
#' 
#' Examples of acceptable inputs:
#' 
#' \code{mutational_effects = list(group_A = mutational_signature_effects(...),
#'                           group_B = mutational_signature_effects(...),
#'                           ...)}
#' 
#' \code{mutational_effects = mutational_signature_effects(...)}                       
#'
#' @param signature_groupings
#' How signatures are grouped in the plot.
#' 1. "auto" shows the superset of the top_n signatures by cancer effect in each
#' group, lumping remaining signatures into "Other".
#'
#' 2. "suggested" supplies signature groupings according to shared etiology
#' as outlined in Cannataro et al. 2022. Also provides colors that convey the etiology
#' of each grouping. If this is used, any provided color palette
#' in "colors" will be overridden.
#'
#' 3. Alternatively, supply a named list where each element contains one or more
#' signature names. The list names will be used to label the groups. If any
#' signatures included in signature_groupings don't exist in the data, they will 
#' not be graphed (but NO warning will be thrown). Remaining signatures will be
#' lumped into "Other".
#' 
#' Example:
#' 
#' \code{signature_groupings = list(clocklike = c('SBS1', 'SBS5'), tobacco = 'SBS4')}
#' 
#' You can also choose to display a manual list of signatures this way:
#' 
#' \code{signature_groupings = list(SBS1 = 'SBS1, SBS5 = 'SBS5')}
#'
#' @param colors
#' Supplies signature group colors.
#' A character vector of colors to use (i.e. a color palette), or
#' "auto". The number of colors provided must be at least the number of
#' signature groupings + 1, for the "Other" category.
#' 
#' @param top_n 
#' How many signatures to include in the plot if signature_groupings = 'auto'.
#' If signature_groupings is not auto, this is not used.

plot_signature_effects = function(mutational_effects = NULL,
                                               signature_groupings = 'auto',
                                               colors = 'auto',
                                               top_n = 5) {
  # Check for proper colors input
  is_color_palette = function(palette_list) {
    valid_colors = colors()
    all_valid = all(sapply(palette_list, function(color) color %in% valid_colors))
    return(all_valid)
  }
  
  if(!(identical(colors, 'auto') |
       is_color_palette(colors)
  )) {
    stop('check that colors is a proper color palette')
  }
  
  # Check that mutational_effects is a proper input, otherwise pass an error.
  # (i.e. a mutational_signature_effects() output, either as a single table or named list)
  if (length(mutational_effects) == 2 &&
      identical(names(mutational_effects),
                c("mutational_sources", "effect_shares"))) {
    # mutational_effects corresponds to expected output (singular table).
    
    # Put it into a list, since the rest of the function uses list syntax.
    # make sample_groupings NA since there's only one group
    mutational_effects = list(mutational_effects)
    sample_groupings = NA
    
  } else if (typeof(mutational_effects) == 'list' &&
             all(sapply(mutational_effects, function(x) length(x) == 2 &&
                        identical(names(x), c("mutational_sources", "effect_shares"))))) {
    # Check that mutational_effects is a named list containing expected outputs.
    # If so, extract names and put into sample_groupings.
    f = function(lst)
      length(lst) == sum(names(lst) != "", na.rm = TRUE)
    if (f(mutational_effects)) {
      # is a named list
      sample_groupings = names(mutational_effects)
    } else {
      stop(
        'check that mutational_effects is a named list or single table of mutational_signature_effects() output(s)'
      )
    }
  } else {
    stop(
      'check that mutational_effects is a named list or single table of mutational_signature_effects() output(s)'
    )
  }
  
  # Create empty list to load data into
  final_df_list = list()
  # Load relevant data from all mutational_signature_effects() outputs into one table for plotting
  for (i in 1:length(mutational_effects)) {
    # extract one mutational_signature_effects() output
    effects_output = mutational_effects[[i]]
    source_share = effects_output$mutational_sources$average_source_shares
    effect_share = effects_output$effect_shares$average_effect_shares
    
    # create separate tables for effect and source share, then merge
    df_weights = data.table(type = "SW",
                            prop = source_share,
                            name = names(source_share))
    df_effects = data.table(type = "CEW",
                            prop = effect_share,
                            name = names(effect_share))
    df = rbind(df_weights, df_effects)
    
    # assign group based on sample_groupings
    df[, sample_group := sample_groupings[[i]]]
    
    # append table to list of SW/CEW data, grouped by sample_groupings
    final_df_list = append(final_df_list, list(df))
  }
  
  # combine all tables
  df_final = do.call(rbind, final_df_list)
  
  # assign signature_groupings to table according to signature_groupings parameter
  if (identical(signature_groupings, 'auto')) {
    # get superset of top_n signatures by CEW in each group
    # first check that top_n is not too big and is not a weird number in general
    if(!(is.numeric(top_n) &
         (0 < top_n) &
         (top_n <= length(mutational_effects[[1]]$mutational_sources$average_source_shares)))) {
      stop('check that top_n is a number between 0 and the number of signatures present in the data')
    }
    
    superset_sigs = list()
    for (i in 1:length(mutational_effects)) {
      mut_effects = mutational_effects[[i]]
      top_sigs = sort(mut_effects$effect_shares$average_effect_shares,
                      decreasing = TRUE)[1:top_n]
      superset_sigs = append(superset_sigs, names(top_sigs))
    }
    
    superset_sigs = unique(unlist(superset_sigs))
    
    # order signatures numerically
    superset_sigs = superset_sigs[order(as.numeric(gsub("\\D", "", superset_sigs)))]
    
    # load into grouping vector
    signature_groupings = list()
    for (sig in superset_sigs) {
      signature_groupings[[sig]] <- sig
    }
  } else if (identical(signature_groupings, 'suggested')) {
    # group signatures into suggested categories as outlined in Cannataro et al.
    signature_groupings = list(
      "Deamination with age, clock-like (1)" = "SBS1",
      "Unknown, clock-like (5)" = "SBS5",
      "APOBEC (2, 13)" = c("SBS2", "SBS13"),
      "Defective homologous recombination (3)" = "SBS3",
      "Tobacco (4, 29)" = c("SBS4", "SBS29"),
      "UV light (7a–d, 38)" = c("SBS7a", "SBS7b", "SBS7c", "SBS7d", "SBS38"),
      "Prior treatment (11, 31, 32, 35)" = c("SBS11", "SBS31", "SBS32", "SBS35"),
      "Mutagenic chemical exposure (22, 24, 42, 88)" = c("SBS22", "SBS24", "SBS42", "SBS88"),
      "Alcohol-associated (16)" = "SBS16"
    )
    
    # set colors to suggested to use custom palette
    colors = 'suggested'
  }
  
  # create group column and assign signature_groupings
  # nested loop syntax is for cases when signature_groupings contains several signatures
  # Ex: "UV light" = SBS7a-d and SBS38
  df_final$group = NA
  
  num_signature_groups = length(signature_groupings)
  for (i in 1:(num_signature_groups)) {
    for (j in 1:(length(signature_groupings[[i]]))) {
      for (k in 1:(nrow(df_final))) {
        if (df_final$name[k] == signature_groupings[[i]][j]) {
          df_final$group[k] = names(signature_groupings[i])
        }
      }
    }
  }
  
  # Use custom palette outlined in Cannataro et al.
  if (identical(colors, 'suggested')) {
    # set all remaining signatures to "Non-actionable..."
    df_final$group = ifelse(is.na(df_final$group),
                            "Non-actionable and unknown signatures",
                            df_final$group)
    
    # Set color palette to match signatures with colors that convey etiology
    palette = c(
      "Deamination with age, clock-like (1)" = "gray40",
      "Unknown, clock-like (5)" = "gray60",
      "APOBEC (2, 13)" = "#7570b3",
      "Defective homologous recombination (3)" = "#e7298a",
      "Tobacco (4, 29)" = "#a6761d",
      "UV light (7a–d, 38)" = "#e6ab02",
      "Prior treatment (11, 31, 32, 35)" = "#1b9e77",
      "Mutagenic chemical exposure (22, 24, 42, 88)" = "#66a61e",
      "Alcohol-associated (16)" = "#d95f02",
      "Non-actionable and unknown signatures" = "black"
    )
  } else if (identical(colors, 'auto')) {
    # set all remaining signatures to 'other'.
    df_final$group = ifelse(is.na(df_final$group),
                            "Other",
                            df_final$group)
  } else {
    # set all remaining signatures to 'other'.
    df_final$group = ifelse(is.na(df_final$group),
                            "Other",
                            df_final$group)
    
    # check that user-provided palette has enough colors
    if(length(colors) < length(unique(df_final$group))) {
      stop(paste0('user-supplied color palette does not have enough colors.
                  (amount needed: ',
                  length(unique(df_final$group)), ')'))
    }
    
    # Use user-provided palette
    palette = colors
  }
  
  # time to plot!
  
  # first, we reorder and rename some values to make the graph pretty.
  # reorder facets because the order is weird
  df_final$facet = factor(df_final$sample_group,
                          levels = sample_groupings)
  # rename values
  df_final$weight_type = ifelse(df_final$type == "SW",
                                "Mutational\n Source Share",
                                "Cancer Effect\n Share")
  # reorder weight bars manually to show Source Share, then Effect Share
  df_final$weight_type = factor(
    df_final$weight_type,
    levels = c("Mutational\n Source Share",
               "Cancer Effect\n Share")
  )
  # reorder signature groups
  if(identical(colors, 'suggested')) {
    # order into same order as in Cannataro et al. 2022
    df_final$group = factor(
      df_final$group,
      levels = c(
        "Non-actionable and unknown signatures",
        "Deamination with age, clock-like (1)",
        "Unknown, clock-like (5)",
        "APOBEC (2, 13)",
        "Defective homologous recombination (3)",
        "Tobacco (4, 29)",
        "UV light (7a–d, 38)",
        "Prior treatment (11, 31, 32, 35)",
        "Mutagenic chemical exposure (22, 24, 42, 88)",
        "Alcohol-associated (16)"
      )
    )
  } else if(exists('superset_sigs')) {
    # order numerically
    df_final$group = factor(
      df_final$group,
      levels = c("Other", superset_sigs)
    )
  } else {
    # force into the same order as original signature_groupings param
    df_final$group = factor(
      df_final$group,
      levels = c("Other", names(signature_groupings))
    )
  }
  
  # run ggplot
  
  if (any(is.na(df_final$facet))) {
    # Single table provided, plot without group argument
    gg = ggplot(data = df_final) +
      geom_bar(
        mapping = aes(
          x = weight_type,
          fill = group,
          weight = prop
        ),
        position = 'fill',
        width = .8,
        color = 'black'
      ) +
      ylab('Weight proportion') +
      xlab('Signature weights and cancer effect weights') +
      theme_classic() + theme(
        axis.title.x = element_text(size = 12),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 12)
      ) +
      scale_y_continuous(n.breaks = 10) + theme(legend.position = 'right')
    
    if(identical(colors, 'auto')) {
      gg = gg + scale_fill_viridis(name = 'Signature', discrete = TRUE) # Using viridis color scale
    } else {
      gg = gg + scale_fill_manual(name = 'Signature', values = rev(palette))
    }
    
  } else {
    # Named list provided, plot with group argument
    gg = ggplot(data = df_final) +
      geom_bar(
        mapping = aes(
          x = weight_type,
          fill = group,
          weight = prop
        ),
        position = 'fill',
        width = .8,
        color = 'black'
      ) +
      ylab('Weight proportion') +
      xlab('Signature weights and cancer effect weights') +
      scale_fill_viridis(discrete = TRUE) + # Using viridis color scale
      theme_classic() + theme(
        axis.title.x = element_text(size = 12),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 12)
      ) +
      scale_y_continuous(n.breaks = 10) + theme(legend.position = 'right') +
      facet_wrap(~ facet, nrow = 1)
    
    if(identical(colors, 'auto')) {
      gg = gg + scale_fill_viridis(name = 'Signature', discrete = TRUE) # Using viridis color scale
    } else {
      gg = gg + scale_fill_manual(name = 'Signature', values = rev(palette))
    }
  }
}

