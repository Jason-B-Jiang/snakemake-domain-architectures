# -----------------------------------------------------------------------------
#
# Plot ortholog length distributions for outgroup and non-outgroup species
#
# Jason Jiang - Created: Apr/12/2023
#               Last edited: Apr/12/2023
#
# Reinke Lab - Microsporidia orthologs
#
# -----------------------------------------------------------------------------

suppressMessages(library(tidyverse))

################################################################################

## Global variabLes (WILL DELETE LATER)

OUTGROUPS <- c('R_allo', 'S_pombe', 'D_disc', 'C_eleg', 'H_sapi', 'D_mela', 'D_reri')
EXCLUDED_SP <- c('P_neur', 'D_roes', 'C_dike', 'N_apis', 'N_bomb', 'H_magn', 'H_tvae', 'D_muel')
LENGTH_LEVELS <- c('<30%', '30% - <60%', '60% - <90%', '≥90%')

################################################################################

main <- function() {
  essential_yeast_genes <- read_lines('../../data/essential_yeast_genes.txt')
  
  length_counts <- read_csv('../../results/single_copy_orthogroups.csv') %>%
    mutate(species_is_outgroup = species %in% OUTGROUPS,
           species_is_outgroup = ifelse(species_is_outgroup, 'Outgroup', 'Microsporidia'),
           essential = ref_ortholog %in% essential_yeast_genes,
           exclude_species = species %in% EXCLUDED_SP,
           length_ratio = ortholog_length / ref_ortholog_length,
           length_bin = get_length_bin(length_ratio)) %>%
    filter(!exclude_species) %>%
    rename(group = species_is_outgroup) %>%
    group_by(group, essential, length_bin) %>%
    summarise(n = n())
  
  # plot stacked bar chart of counts of orthologs in each length bin for
  # Microsporidia and outgroups, stacked by essentiality of orthologs
  ggplot(length_counts, aes(x = factor(length_bin, levels = LENGTH_LEVELS),
                            y = n, fill = essential)) +
    geom_bar(position = 'dodge', stat = 'identity', color = 'black') +
    geom_text(aes(label = n), position = position_dodge(width = 0.9), vjust = -0.5) +
    labs(x = '% Length to Saccharomyces cerevisiae ortholog', y = 'Count', fill = 'Essential') +
    scale_fill_manual('Ortholog essential in yeast?',
                      values = c('#F8766D', '#619CFF')) + # set fill colors
    facet_wrap(~group, ncol = 2) + # add facet wrap by 'group'
    theme_minimal() +
    theme(
      axis.text = element_text(size = 14, color = 'black'), # set color and font size of axis text
      axis.title = element_text(size = 14, color = 'black'), # set color and font size of axis titles
      legend.text = element_text(size = 14, color = 'black'), # set color and font size of legend text
      strip.text = element_text(size = 14, color = 'black'), # set color and font size of facet wrap titles
      legend.justification = 'center',
      legend.position = 'bottom',
      legend.title = element_text(size = 14),
      title = element_text(color = 'black', size = 14)
    )
}

################################################################################

## Helper functions

get_length_bin <- Vectorize(function(length_ratio) {
  if (length_ratio < 0.30) {
    return('<30%')
  } else if (length_ratio >= 0.30 & length_ratio < 0.60) {
    return('30% - <60%')
  } else if (length_ratio >= 0.60 & length_ratio < 0.90) {
    return('60% - <90%')
  } else {
    return('≥90%')
  }
})

################################################################################

main()