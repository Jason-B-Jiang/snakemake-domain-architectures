# -----------------------------------------------------------------------------
#
# Do pairwise alignment of ortholog domain architectures within orthogroups
#
# Jason Jiang - Created: Feb/28/2023
#               Last edited: Feb/28/2023
#
# Reinke Lab - Microsporidia orthologs
#
# -----------------------------------------------------------------------------

library(tidyverse)
library(NameNeedle)

################################################################################

## Constants

# Column names for cath-resolve-hits output files
CRH_HEADER = c('query_id', 'match_id', 'score', 'boundaries', 'resolved',
               'cond_evalue', 'indp_evalue')

# Parameters for Needleman-Wunsch alignments of domain architectures
NW_PARAMS = list(MATCH = 0, MISMATCH = -3, GAP = -10, GAPCHAR = '*')

# List of abbreviated outgroup species names
# TODO - integrate with Snakemake config file
OUTGROUPS <- c('C_eleg', 'D_disc', 'D_reri', 'D_mela', 'H_sapi', 'S_pomb')

################################################################################

main <- function() {
  # args <- commandArgs(trailingOnly = T)
  domain_archs_dir <- "../../results/domain_architectures"
  pfam_clans <- make_pfam_clan_hashtable("../../data/pfam/Pfam-A-clans.tsv")
  out <- '../../results/aligned_domain_architectures.csv'
  
  merged_domain_archs <- merge_domain_archs(domain_archs_dir, pfam_clans)
  
  write_csv(aligned_domain_archs, out)
}

################################################################################

## Helper functions

merge_domain_archs <- function(domain_archs_dir, pfam_clans) {
  # ---------------------------------------------------------------------------
  # Docstring goes here.
  # ---------------------------------------------------------------------------
  domain_archs_df <- data.frame(matrix(ncol = 6, nrow = 0))
  colnames(domain_archs_df) <- c('orthogroup', 'species', 'query_id', 'match_id',
                                 'resolved', 'match_clans')
  
  for (orthogroup_folder in list.files(domain_archs_dir, full.names = TRUE)) {
    orthogroup <- basename(orthogroup_folder)
    
    for (crh_file in list.files(orthogroup_folder, full.names = TRUE)) {
      crh_df <- read_delim(crh_file, delim = ' ', comment = '#',
                           col_names = CRH_HEADER, show_col_types = F)
      
      if (nrow(crh_df) > 0) {
        domain_archs_df <- rbind(domain_archs_df,
                                 parse_crh_output(crh_df, pfam_clans, orthogroup))
      }
    }
  }
  
  return(domain_archs_df)
}


make_pfam_clan_hashtable <- function(pfam_clans_file) {
  # ---------------------------------------------------------------------------
  # Docstring goes here.
  # ---------------------------------------------------------------------------
  pfam_clans_df <- read_tsv(pfam_clans_file, show_col_types = F)
  
  # create hashtables mapping both pfam families back to their clans
  pfam_clans <- new.env()
  Map(function(fam, clan) {pfam_clans[[fam]] <- ifelse(is.na(clan), fam, clan)},
      pfam_clans_df$Family_ID, pfam_clans_df$Clan_name)
  
  return(pfam_clans)
}


parse_crh_output <- function(crh_df, pfam_clans, orthogroup) {
  # ---------------------------------------------------------------------------
  # Read in cath-resolve-hits output file for an orthogroup's ortholog domain
  # architectures + domain boundaries
  #
  # Args:
  #   resolved_hits: filepath to cath-resolve-hits output file
  #
  # ---------------------------------------------------------------------------
  return(
    crh_df %>%
      select(query_id, match_id, resolved) %>%
      mutate(orthogroup = orthogroup,
             species = NA) %>%
      group_by(query_id) %>%
      mutate(match_id = str_c(match_id, collapse = '; '),
             resolved = str_c(resolved, collapse = '; ')) %>%
      ungroup() %>%
      rowwise() %>%
      mutate(match_clans = get_domain_clans(match_id, pfam_clans)) %>%
      ungroup() %>%
      group_by(query_id) %>%
      distinct(.keep_all = T)
  )
}


get_domain_clans <- function(domains, pfam_clans) {
  # ---------------------------------------------------------------------------
  # Docstring goes here
  # ---------------------------------------------------------------------------
  domains <- str_split(domains, '; ')[[1]]
  
  return(str_c(
    # if domain isn't in pfam_clans hashtable, just keep the domain as is
    sapply(domains, function(x) {ifelse(!is.null(pfam_clans[[x]]),
                                        pfam_clans[[x]],
                                        x)}),
    collapse = '; '))
}

################################################################################

main()