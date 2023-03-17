# -----------------------------------------------------------------------------
#
# Do pairwise alignment of ortholog domain architectures within orthogroups
#
# Jason Jiang - Created: Feb/28/2023
#               Last edited: Mar/17/2023
#
# Reinke Lab - Microsporidia orthologs
#
# -----------------------------------------------------------------------------

suppressMessages(library(tidyverse))
suppressMessages(library(NameNeedle))

################################################################################

## Constants

# Column names for cath-resolve-hits output files
CRH_HEADER = c('ortholog', 'domain_arch', 'score', 'boundaries', 'domain_bounds',
               'cond_evalue', 'indp_evalue')

# Parameters for Needleman-Wunsch alignments of domain architectures
NW_PARAMS = list(MATCH = 0, MISMATCH = -3, GAP = -10, GAPCHAR = '*')

################################################################################

main <- function() {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  args <- commandArgs(trailingOnly = T)
  domain_archs_dir <- args[1]
  orthogroups <- read_csv(args[2], show_col_types = FALSE)
  pfam_clans <- make_pfam_clan_hashtable(args[3])
  out <- args[4]
  
  # domain_archs_dir <- "../../results/domain_architectures"
  # orthogroups <- read_csv("../../results/single_copy_orthogroups.csv",
  #                         show_col_types = FALSE)
  # pfam_clans <- make_pfam_clan_hashtable("../../data/pfam/Pfam-A-clans.tsv")
  # out <- '../../results/aligned_domain_architectures.csv'

  aligned_domain_archs <- merge_domain_archs(domain_archs_dir,
                                             orthogroups,
                                             pfam_clans) %>%
    rowwise() %>%
    mutate(aligned_domain_archs = align_domain_archs(domain_arch_clans,
                                                     ref_domain_arch_clans),
           lost_doms = get_domain_architecture_diffs(aligned_domain_archs)[['loss']],
           gained_doms = get_domain_architecture_diffs(aligned_domain_archs)[['gain']],
           swapped_doms = get_domain_architecture_diffs(aligned_domain_archs)[['swap']])
  
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

merge_domain_archs <- function(domain_archs_dir, orthogroups, pfam_clans) {
  # ---------------------------------------------------------------------------
  # Docstring goes here.
  # ---------------------------------------------------------------------------
  domain_archs_df <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(domain_archs_df) <- c('orthogroup', 'ortholog', 'domain_arch',
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
                                                ortholog_name))
    }
  }
  
  species_domain_archs <- domain_archs_df %>%
    filter(!(ortholog %in% orthogroups$ref_ortholog)) %>%
    right_join(select(orthogroups, species, is_outgroup, ortholog,
                      ortholog_length), by = 'ortholog')
  
  ref_domain_archs <- domain_archs_df %>%
    filter(ortholog %in% orthogroups$ref_ortholog) %>%
    rename(ref_ortholog = ortholog, ref_domain_arch = domain_arch,
           ref_domain_bounds = domain_bounds,
           ref_domain_arch_clans = domain_arch_clans) %>%
    right_join(select(orthogroups, ref_species, ref_ortholog,
                      ref_ortholog_length), by = 'ref_ortholog') %>%
    distinct(.keep_all = TRUE)
  
  return(
    full_join(species_domain_archs,
              ref_domain_archs,
              by = 'orthogroup') %>%
      select(orthogroup, species, is_outgroup, ref_species,
             ortholog, ref_ortholog, ortholog_length, ref_ortholog_length,
             domain_arch, ref_domain_arch, domain_arch_clans,
             ref_domain_arch_clans, domain_bounds, ref_domain_bounds)
  )
}


parse_crh_output <- function(crh_df, pfam_clans, orthogroup, ortholog_name) {
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
      mutate(domain_arch_clans = get_domain_clans(domain_arch, pfam_clans)) %>%
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


get_clan_letters <- function(clans) {
  # ---------------------------------------------------------------------------
  # Get unique letter mappings for a vector of clans
  # ---------------------------------------------------------------------------
  clans_to_letters <- lapply(as.list(1 : length(clans)), function(i) {LETTERS[i]})
  names(clans_to_letters) <- clans
  
  return(clans_to_letters)
}


format_alignment_str <- function(aligned_DA, DA) {
  # ---------------------------------------------------------------------------
  # Replace letters from a Needleman-Wunsch domain architecture alignment with
  # the Pfam clans from the original domain architecture.
  # ---------------------------------------------------------------------------
  aligned_DA <- str_split(aligned_DA, '')[[1]]
  
  i = 1
  for (j in 1 : length(aligned_DA)) {
    if (aligned_DA[j] != '*') {
      aligned_DA[j] <- DA[i]
      i <- i + 1
    }
  }
  
  return(str_c(aligned_DA, collapse = ' -> '))
}


get_domain_architecture_diffs <- function(aligned_domain_archs) {
  # ---------------------------------------------------------------------------
  # Get domain architectural changes (loss, gain, swap) in an ortholog relative
  # to a reference ortholog.
  #
  # Return the losses (loss, gain, swap) as a named list of vectors.
  # ---------------------------------------------------------------------------
  aligned_domain_archs <- str_split(aligned_domain_archs, '; ')[[1]]
  domain_arch <- str_split(aligned_domain_archs[1], ' -> ')[[1]]
  ref_domain_arch <- str_split(aligned_domain_archs[2], ' -> ')[[1]]
  
  loss <- which(domain_arch == '*')
  gain <- which(ref_domain_arch == '*')
  swap <- which(domain_arch != '*' & ref_domain_arch != '*' & domain_arch != ref_domain_arch)
  
  # replace indices of lost/swapped/gained domains with the actual domains
  loss <- ifelse(length(loss) > 0,
                 str_c(ref_domain_arch[loss], collapse = '; '),
                 NA)
  
  gain <- ifelse(length(gain) > 0,
                 str_c(domain_arch[gain], collapse = '; '),
                 NA)
  
  swap <- ifelse(length(swap) > 0,
                 str_c(mapply(function(x, y) {str_c(x, y, sep = ' ~> ')},
                        ref_domain_arch[swap],
                        domain_arch[swap]), sep = '; '),
                 NA)
  
  return(list('loss' = loss, 'gain' = gain, 'swap' = swap))
}

################################################################################

main()