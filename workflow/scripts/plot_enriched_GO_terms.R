# -----------------------------------------------------------------------------
#
# Do pairwise alignment of ortholog domain architectures within orthogroups
#
# Jason Jiang - Created: Feb/28/2023
#               Last edited: Mar/27/2023
#
# Reinke Lab - Microsporidia orthologs
#
# -----------------------------------------------------------------------------

suppressMessages(library(tidyverse))

################################################################################

main <- function() {
  go_bio <- read_tsv('../../results/go_enrichment/go_bio.txt') %>%
    filter(go_term != 'Unclassified (UNCLASSIFIED)',
           fold_enrich == '+') %>%
    mutate(go_id = str_extract(go_term, 'GO:\\d+')) %>%
    select(go_id, go_term, everything())
  
  go_molec <- read_tsv('../../results/go_enrichment/go_molec.txt') %>%
    filter(go_term != 'Unclassified (UNCLASSIFIED)',
           fold_enrich == '+') %>%
    mutate(go_id = str_extract(go_term, 'GO:\\d+')) %>%
    select(go_id, go_term, everything())
  
  go_comp <- read_tsv('../../results/go_enrichment/go_comp.txt') %>%
    filter(go_term != 'Unclassified (UNCLASSIFIED)',
           fold_enrich == '+') %>%
    mutate(go_id = str_extract(go_term, 'GO:\\d+')) %>%
    select(go_id, go_term, everything())
  
  # Write GO IDs and FDRs as text files to copy into revigo
  write_lines(str_c(go_bio$go_id, go_bio$fdr, sep = '  '), '../../results/go_enrichment/revigo_bio.txt')
  write_lines(str_c(go_molec$go_id, go_molec$fdr, sep = '  '), '../../results/go_enrichment/revigo_molec.txt')
  write_lines(str_c(go_comp$go_id, go_comp$fdr, sep = '  '), '../../results/go_enrichment/revigo_comp.txt')
  
  # Load in trimmed GO terms from revigo
  revigo_bio <- read_tsv('../../results/go_enrichment/Revigo_BP_OnScreenTable.tsv') %>%
    filter(Representative == 'null')
  
  revigo_molec <- read_tsv('../../results/go_enrichment/Revigo_MF_OnScreenTable.tsv') %>%
    filter(Representative == 'null')
  
  revigo_comp <- read_tsv('../../results/go_enrichment/Revigo_CC_OnScreenTable.tsv') %>%
    filter(Representative == 'null')
  
  # filter GO terms to only trimmed GO terms from revigo, and sort in increasing
  # order of FDR, limiting to top 20 GO terms only
  go_bio <- go_bio %>%
    filter(go_id %in% revigo_bio$TermID) %>%
    arrange(fdr) %>%
    top_n(-20) %>%  # choose top 20 significant terms
    arrange(-enrich_direction) %>%  # sort in decreasing order of fold enrichment
    rename('Genes annotated' = n_observed, 'Adjusted p-value (BH)' = fdr,
           'Fold enrichment' = enrich_direction, 'GO term' = go_term) %>%
    mutate(`-log₂ ( Adjusted p-value (BH) )` = -log(`Adjusted p-value (BH)`))
  
  go_molec <- go_molec %>%
    filter(go_id %in% revigo_molec$TermID) %>%
    arrange(fdr) %>%
    top_n(-20) %>%  # choose top 20 significant terms
    arrange(-enrich_direction) %>%  # sort in decreasing order of fold enrichment
    rename('Genes annotated' = n_observed, 'Adjusted p-value (BH)' = fdr,
           'Fold enrichment' = enrich_direction, 'GO term' = go_term) %>%
    mutate(`-log₂ ( Adjusted p-value (BH) )` = -log(`Adjusted p-value (BH)`))
  
  go_comp <- go_comp %>%
    filter(go_id %in% revigo_comp$TermID) %>%
    arrange(fdr) %>%
    top_n(-20) %>%  # choose top 20 significant terms
    arrange(-enrich_direction) %>%  # sort in decreasing order of fold enrichment
    rename('Genes annotated' = n_observed, 'Adjusted p-value (BH)' = fdr,
           'Fold enrichment' = enrich_direction, 'GO term' = go_term) %>%
    mutate(`-log₂ ( Adjusted p-value (BH) )` = -log(`Adjusted p-value (BH)`))
  
  # plot the enriched GO terms
  fdr_range <- range(rbind(go_bio, go_comp, go_molec)[['-log₂ ( Adjusted p-value (BH) )']])
  count_range <- range(rbind(go_bio, go_comp, go_molec)[['Genes annotated']])
  
  go_bio_plot <- ggplot(go_bio, aes(x = `Fold enrichment`,
                                    y = rev(factor(`GO term`, go_bio$`GO term`)),
                                    color = `-log₂ ( Adjusted p-value (BH) )`,
                                    size = `Genes annotated`)) +
    geom_point() +
    scale_color_gradient(low = 'red', high = 'blue', limits = fdr_range) +
    scale_size_continuous(range = c(2, 10), limits = count_range) +
    theme_bw() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(color = 'black', size = 14),
          axis.title.x = element_text(color = 'black', size = 18),
          axis.text.x = element_text(color = 'black', size = 14),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 18))
  
  go_molec_plot <- ggplot(go_molec, aes(x = `Fold enrichment`,
                                    y = rev(factor(`GO term`, go_molec$`GO term`)),
                                    color = `Adjusted p-value (BH)`,
                                    size = `Genes annotated`)) +
    geom_point() +
    scale_color_gradient(low = 'red', high = 'blue', limits = fdr_range) +
    scale_size_continuous(range = c(2, 10), limits = count_range) +
    theme_bw() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(color = 'black', size = 14),
          axis.title.x = element_text(color = 'black', size = 18),
          axis.text.x = element_text(color = 'black', size = 14),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 18))
  
  go_comp_plot <- ggplot(go_comp, aes(x = `Fold enrichment`,
                                    y = rev(factor(`GO term`, go_comp$`GO term`)),
                                    color = `Adjusted p-value (BH)`,
                                    size = `Genes annotated`)) +
    geom_point() +
    scale_color_gradient(low = 'red', high = 'blue', limits = fdr_range) +
    scale_size_continuous(range = c(2, 10), limits = count_range) +
    theme_bw() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(color = 'black', size = 14),
          axis.title.x = element_text(color = 'black', size = 18),
          axis.text.x = element_text(color = 'black', size = 14),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 18))
}

################################################################################

main()