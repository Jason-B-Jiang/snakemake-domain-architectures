# -----------------------------------------------------------------------------
#
# Do pairwise alignment of ortholog domain architectures within orthogroups
#
# Jason Jiang - Created: Feb/28/2023
#               Last edited: Mar/07/2023
#
# Reinke Lab - Microsporidia orthologs
#
# -----------------------------------------------------------------------------

library(tidyverse)
library(NameNeedle)

################################################################################

## Constants

# Column names for cath-resolve-hits output files
CRH_HEADER = c('ortholog', 'domain_arch', 'score', 'boundaries', 'domain_bounds',
               'cond_evalue', 'indp_evalue')

# Parameters for Needleman-Wunsch alignments of domain architectures
NW_PARAMS = list(MATCH = 0, MISMATCH = -3, GAP = -10, GAPCHAR = '*')

################################################################################

main <- function() {
  args <- commandArgs(trailingOnly = T)
  
  # domain_archs_dir <- args[1]
  # pfam_clans <- make_pfam_clan_hashtable(args[2])
  # ortho_to_species <- make_ortholog_to_species_hashtable(
  #   read_tsv(args[3], show_col_types = FALSE)
  #   )
  # out <- args[4]
  # reference_species <- args[5]
  # outgroups <- args[6 : length(args)]
  
  domain_archs_dir <- "../../results/domain_architectures"
  pfam_clans <- make_pfam_clan_hashtable("../../data/pfam/Pfam-A-clans.tsv")
  ortho_to_species <- make_ortholog_to_species_hashtable(
    read_tsv('../../results/OrthoFinder/Results_OrthoFinder/Orthogroups/Orthogroups.tsv')
  )
  out <- '../../results/aligned_domain_architectures.csv'
  reference_species <- 'E_brev'
  outgroups <- c('E_cuni')

  aligned_domain_archs <- merge_domain_archs(domain_archs_dir,
                                             pfam_clans,
                                             ortho_to_species) %>%
    get_pairwise_domain_archs(reference_species) %>%
    mutate(aligned_domain_archs = align_domain_archs(domain_arch_clans,
                                                     ref_domain_arch_clans),
           lost_doms = NA,
           swapped_doms = NA,
           gained_doms = NA)
  
  write_csv(aligned_domain_archs, out)
}

################################################################################

## Helper functions

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


make_ortholog_to_species_hashtable <- function(orthogroups) {
  # ---------------------------------------------------------------------------
  # Docstring goes here.
  # ---------------------------------------------------------------------------
  ortho_to_species <- new.env()
  
  for (species in colnames(orthogroups)[2 : ncol(orthogroups)]) {
    orthologs <-
      unlist(flatten(lapply(orthogroups[[species]][!is.na(orthogroups[[species]])], function(x) {str_split(x, ', ')[[1]]})))
    
    for (ortholog in orthologs) {
      ortho_to_species[[ortholog]] <- species
    }
  }
  
  return(ortho_to_species)
}


merge_domain_archs <- function(domain_archs_dir, pfam_clans, ortho_to_species) {
  # ---------------------------------------------------------------------------
  # Docstring goes here.
  # ---------------------------------------------------------------------------
  domain_archs_df <- data.frame(matrix(ncol = 6, nrow = 0))
  colnames(domain_archs_df) <- c('orthogroup', 'species', 'ortholog', 'domain_arch',
                                 'domain_bounds', 'domain_arch_clans')
  
  for (orthogroup_folder in list.files(domain_archs_dir, full.names = TRUE)) {
    orthogroup <- basename(orthogroup_folder)
    
    for (crh_file in list.files(orthogroup_folder, full.names = TRUE)) {
      crh_df <- read_delim(crh_file, delim = ' ', comment = '#',
                           col_names = CRH_HEADER, show_col_types = F)
      
      ortholog_name <- basename(crh_file)
      
      domain_archs_df <- rbind(domain_archs_df,
                               parse_crh_output(crh_df,
                                                pfam_clans,
                                                orthogroup,
                                                ortho_to_species,
                                                ortholog_name))
    }
  }
  
  # arrange columns in desired order before returning dataframe, and filter out
  # orthologs in orthogroups that are not single-copy
  return(
    domain_archs_df %>%
      group_by(orthogroup, species) %>%
      # filter out non single-copy orthologs from orthogroups
      filter(n() < 2) %>%
      ungroup() %>%
      select(orthogroup, species, ortholog, domain_arch, domain_arch_clans,
             domain_bounds)
  )
}


parse_crh_output <- function(crh_df, pfam_clans, orthogroup, ortho_to_species,
                             ortholog_name) {
  # ---------------------------------------------------------------------------
  # Read in cath-resolve-hits output file for an orthogroup's ortholog domain
  # architectures + domain boundaries
  #
  # Args:
  #   resolved_hits: filepath to cath-resolve-hits output file
  #
  # ---------------------------------------------------------------------------
  # Empty cath-resolve-hits file for ortholog, so no domain assignments were
  # made
  if (nrow(crh_df) < 1) {
    return(data.frame(orthogroup = orthogroup,
                      species = ortho_to_species[[ortholog_name]],
                      ortholog = ortholog_name,
                      domain_arch = NA,
                      domain_bounds = NA,
                      domain_arch_clans = NA))
  }
  
  return(
    crh_df %>%
      select(ortholog, domain_arch, domain_bounds) %>%
      mutate(orthogroup = orthogroup) %>%
      group_by(ortholog) %>%
      mutate(domain_arch = str_c(domain_arch, collapse = '; '),
             domain_bounds = str_c(domain_bounds, collapse = '; ')) %>%
      ungroup() %>%
      rowwise() %>%
      mutate(domain_arch_clans = get_domain_clans(domain_arch, pfam_clans),
             species = ortho_to_species[[ortholog]]) %>%
      ungroup() %>%
      group_by(ortholog) %>%
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


get_pairwise_domain_archs <- function(merged_domain_archs, reference_species) {
  # ---------------------------------------------------------------------------
  # Docstring goes here
  # ---------------------------------------------------------------------------
  # domain architectures for non-reference species
  non_ref <- merged_domain_archs %>%
    filter(species != reference_species)
  
  # domain architectures for reference species
  ref <- merged_domain_archs %>%
    filter(species == reference_species) %>%
    rename(ref_species = species,
           ref_ortholog = ortholog,
           ref_domain_arch = domain_arch,
           ref_domain_arch_clans = domain_arch_clans,
           ref_domain_bounds = domain_bounds)
  
  # do a full join between the non-reference and reference species domain
  # architecture dataframes, to get pairs of reference + non-reference species
  # orthologs
  #
  # NOTE: we are guaranteed to have a single-copy ortholog from the reference
  # species in each orthogroup, so this will be a many-to-one join
  return(full_join(non_ref, ref, by = 'orthogroup'))
}


align_domain_archs <- function(DA_1, DA_2) {
  # ---------------------------------------------------------------------------
  # Align domain architectures between two proteins, using the clans for their
  # domain architectures
  # ---------------------------------------------------------------------------
  DA_1 <- str_split(DA_1, '; ')[[1]]
  DA_2 <- str_split(DA_2, '; ')[[1]]
  
  if (all(is.na(DA_1), is.na(DA_2))) {
    # no domain assignments for either ortholog
    return(NA)
    
    # only domain assignments for one ortholog
  } else if (all(is.na(DA_1))) {
    return(
      str_c(str_c(rep('*', length(DA_2)), collapse = ' -> '),
            str_c(DA_2, collapse = ' -> '),
            sep = '; ')
    )
    
  } else if (all(is.na(DA_2))) {
    return(
      str_c(str_c(DA_1, collapse = ' -> '),
            str_c(rep('*', length(DA_1)), collapse = ' -> '),
            sep = '; ')
    )
  }
  
  # domain architectures available for both orthologs, so align normally
  clans_to_letters <- get_clan_letters(unique(c(DA_1, DA_2)))
  
  # create string representations of each domain architecture using the clan to
  # letter mappings, for alignment
  DA_1_str = str_c(sapply(DA_1, function(x) {clans_to_letters[[x]]}),
                   collapse = '')
  
  DA_2_str = str_c(sapply(DA_2, function(x) {clans_to_letters[[x]]}),
                   collapse = '')
  
  # align domain architectures using their string representations
  aligned_DA <- needles(DA_1_str, DA_2_str, params = NW_PARAMS)
  
  # return formatted string of domain architecture alignments, with original
  # domains replacing the letter representations in the alignments
  return(
    str_c(
      format_alignment_str(aligned_DA$align1, DA_1),
      format_alignment_str(aligned_DA$align2, DA_2),
      sep = '; '
    )
  )
}


format_alignment_str <- function(aligned_DA, DA) {
  return(NA)
}


get_clan_letters <- function(clans) {
  return()
}

################################################################################

main()