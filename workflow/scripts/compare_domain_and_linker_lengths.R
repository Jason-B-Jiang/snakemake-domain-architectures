# -----------------------------------------------------------------------------
#
# Compare domain versus linker amino acid conservation
#
# Jason Jiang - Created: Apr/13/2023
#               Last edited: Apr/13/2023
#
# Reinke Lab - Microsporidia orthologs
#
# -----------------------------------------------------------------------------

suppressMessages(library(tidyverse))
suppressMessages(library(ggsignif))

################################################################################

main <- function() {
  ref_essential_genes <- read_lines('../../data/essential_yeast_genes.txt')
  
  domain_linker_length_ratios <- 
    read_csv('../../results/aligned_ortholog_domain_architectures.csv') %>%
    filter(!exclude_species, !species_is_outgroup, !is.na(aligned_domain_archs),
           ortholog_length < ref_ortholog_length, species != 'R_allo') %>%
    mutate(total_domain_length = get_total_domain_length(domain_lengths),
           ref_total_domain_length = get_total_domain_length(ref_domain_lengths),
           linker_length = ortholog_length - total_domain_length,
           ref_linker_length = ref_ortholog_length - ref_total_domain_length,
           'Domains' = ref_total_domain_length / total_domain_length,
           'Linkers' = ref_linker_length / linker_length,
           essential_in_ref = ref_ortholog %in% ref_essential_genes) %>%
    filter(!is.infinite(`Domains`), !is.infinite(`Linkers`),
           !is.nan(`Domains`), !is.nan(`Linkers`))
  
  # p-value calculations
  domain_vs_linker_essential <- wilcox.test(  # p = 0.0033, bonf corrected
    filter(domain_linker_length_ratios, essential_in_ref)[['Domains']],
    filter(domain_linker_length_ratios, essential_in_ref)[['Linkers']],
    paired = TRUE
  )$p.value
  
  domain_vs_linker_nonessential <- wilcox.test(  # 6.49e-17, bonf corrected
    filter(domain_linker_length_ratios, !essential_in_ref)[['Domains']],
    filter(domain_linker_length_ratios, !essential_in_ref)[['Linkers']],
    paired = TRUE
  )$p.value
  
  # reformat dataframe for plotting
  domain_linker_length_ratios <- domain_linker_length_ratios %>%
    pivot_longer(cols = c('Domains', 'Linkers'),
                 names_to = 'type',
                 values_to = 'length_ratio') %>%
    mutate(sqrt_length_ratio = sqrt(length_ratio))
  
  ggplot(data = domain_linker_length_ratios,
         aes(x = type, y = sqrt_length_ratio, fill = essential_in_ref)) +
    geom_violin() +
    scale_fill_manual('Ortholog essential in yeast?',
                      values = c('#F8766D', '#619CFF')) +
    labs(title = str_c('n = ', nrow(domain_linker_length_ratios) / 2,
                       ' microsporidia-yeast single-copy ortholog pairs\n')) +
    labs(y = 'Length in yeast ortholog / length in microsporidia ortholog (sqrt)') +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size = 14, color = 'black', angle = 45,
                                     vjust = 1, hjust = 1),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 14, color = 'black'),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          legend.position = 'bottom',
          legend.justification = 'center',
          title = element_text(color = 'black', size = 14))
}

################################################################################

## Helper functions

get_total_domain_length <- Vectorize(function(domain_lengths) {
  if (is.na(domain_lengths)) {
    return(0)
  }
  
  return(sum(as.integer(str_split(domain_lengths, '; ')[[1]])))
})

################################################################################

main()