# -----------------------------------------------------------------------------
#
# Plot frequencies of domain architectural change events in orthologs across
# clades to reference species (yeast)
#
# Jason Jiang - Created: 2022/09/19
#               Last edited: 2023/04/03
#
# Reinke Lab - Microsporidia orthologs
#
# -----------------------------------------------------------------------------

library(tidyverse)

################################################################################

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  essential_ref_genes <- read_lines('../../data/essential_yeast_genes.txt')
  domain_arch_changes <- read_csv('../../results/aligned_ortholog_domain_architectures.csv') %>%
    filter(!exclude_species,
           is.na(short_enough_for_domain_loss) | short_enough_for_domain_loss,
           !is.na(aligned_domain_archs)) %>%
    mutate(ref_ortholog_essential = ref_ortholog %in% essential_ref_genes,
           DA_conservation = get_DA_conservation(lost_doms, gained_doms, swapped_doms),
           species_is_outgroup = ifelse(species_is_outgroup, 'Outgroup', 'Microsporidia')) %>%
    separate_rows(DA_conservation, sep = '; ') %>%
    select(species, ortholog, species_is_outgroup, ref_ortholog_essential, DA_conservation)
}

################################################################################

## Helper functions

make_clades_hash <- function(species_clades) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  clades <- new.env()
  map2_chr(species_clades$species,
           species_clades$clade,
           function(x, y) {clades[[x]] <- y})
  
  return(clades)
}


get_DA_conservation <- Vectorize(function(lost_doms, gained_doms, swapped_doms) {
  if (is.na(lost_doms) & is.na(gained_doms) & is.na(swapped_doms)) {
    return("Conserved domain architecture")
  }
  
  DA_changes <- c()
  
  if (!is.na(lost_doms)) {
    DA_changes <- c(DA_changes, 'Lost domain(s)')
  }
  
  if (!is.na(gained_doms)) {
    DA_changes <- c(DA_changes, 'Gained domain(s)')
  }
  
  if (!is.na(swapped_doms)) {
    DA_changes <- c(DA_changes, 'Swapped domain(s)')
  }
  
  return(str_c(DA_changes, collapse = '; '))
}, vectorize.args = c('lost_doms', 'gained_doms', 'swapped_doms'))


plot_DA_change_rates <- function(domain_arch_changes) {
  ggplot(data = domain_arch_changes, aes(x = species_is_outgroup,
                                         fill = DA_conservation)) +
    geom_bar(position = 'dodge')
}

################################################################################

main()