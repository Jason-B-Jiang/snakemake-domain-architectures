# -----------------------------------------------------------------------------
#
# Get single-copy orthogroups with respect to a reference species
#
# Jason Jiang - Created: Feb/01/2023
#               Last edited: Feb/01/2023
#
# Reinke Lab - Microsporidia orthologs
#
# -----------------------------------------------------------------------------

suppressMessages(library(tidyverse))

################################################################################

main <- function() {
  # ---------------------------------------------------------------------------
  # Command line arguments:
  # $1 = filepath to directory with OrthoFinder output
  # $2 = name of reference species to find single-copy orthogroups for
  # $3 = filepath to rds of ortholog lengths in each orthogroup
  # $4 = filepath to save single-copy orthogroups to
  # $5 = filepath to save NAMES of single-copy orthogroups to
  #      (this will help later when choosing orthogroups to run hmmscan on)
  # $6 - $n = names of outgroup species
  # ---------------------------------------------------------------------------
  args <- commandArgs(trailingOnly = TRUE)
  orthogroups <- read_tsv(str_c(args[1], '/Results_OrthoFinder/Orthogroups/Orthogroups.tsv'),
                          show_col_types = FALSE)
  ref_species <- args[2]
  ortho_lengths <- read_rds(args[3])
  orthogroup_out <- args[4]
  orthogroup_names_out <- args[5]
  outgroups <- args[6 : length(args)]
  
  rows_to_include <- c()
  
  # keep rows (i.e: orthogroups) that have a single-copy ortholog from the
  # reference species, and a single-copy ortholog from any other species
  # (yeah, couldn't figure out a nice way to do this with dplyr...)
  for (i in 1 : nrow(orthogroups)) {
    ref_ortholog <- orthogroups[[ref_species]][i]
    
    # orthogroup has single-copy ortholog for reference species
    if (!is.na(ref_ortholog) & !str_detect(ref_ortholog, ',')) {
      
      # check for single-copy ortholog from other species
      for (col in setdiff(colnames(orthogroups), c(ref_species, 'Orthogroup'))) {
        other_ortholog <- orthogroups[[col]][i]
        
        # at least one single-copy ortholog exists from another species in this
        # orthogroup, so consider this orthogroup
        if (!is.na(other_ortholog) & !str_detect(other_ortholog, ',')) {
          rows_to_include <- c(rows_to_include, i)
          break
        }
      }
    }
  }
  
  # filter to orthogroups meeting the criterion mentioned above
  orthogroups <- orthogroups[rows_to_include,]
  
  # format orthogroups for saving as a dataframe, in which row is a single
  # copy ortholog pair between the reference species and another species,
  # long with the lengths for each ortholog
  orthogroups <- format_orthogroups(orthogroups,
                                    ref_species,
                                    ortho_lengths,
                                    outgroups)
  
  # write orthogroups as csv file to specified out filepath
  write_csv(orthogroups, orthogroup_out)
  
  # write text file of orthogroups fitting this criterion to resources folder
  # (will need when picking which orthogroups to assign domains to with hmmscan)
  write_lines(unique(orthogroups$orthogroup), orthogroup_names_out)
}

################################################################################

## Helper functions

format_orthogroups <- function(orthogroups, ref_species, ortho_lengths, outgroups) {
  # ----------------------------------------------------------------------------
  # Docstring goes here.
  # ----------------------------------------------------------------------------
  species_orthos <- select(orthogroups, -ref_species) %>%
    rename(orthogroup = Orthogroup)
  
  ref_orthos <- select(orthogroups, Orthogroup, ref_species) %>%
    rename(orthogroup = Orthogroup) %>%
    mutate(ref_species = ref_species)
  
  colnames(ref_orthos)[2] = 'ref_ortholog'
  
  return(
    species_orthos %>%
      pivot_longer(cols = colnames(species_orthos)[2 : ncol(species_orthos)],
                   names_to = "species",
                   values_to = "ortholog") %>%
      full_join(ref_orthos, by = 'orthogroup') %>%
      # remove entries with missing ortholog for species, or non-single copy ortholog
      filter(!is.na(ortholog), !str_detect(ortholog, ',')) %>%
      mutate(is_outgroup = species %in% outgroups) %>%
      rowwise() %>%
      mutate(ortholog_length = ortho_lengths[[species]][[ortholog]],
             ref_ortholog_length = ortho_lengths[[ref_species]][[ref_ortholog]]) %>%
      # re-order columns to look nicer
      select(orthogroup, species, is_outgroup, ref_species, ortholog,
             ref_ortholog, ortholog_length, ref_ortholog_length)
  )
}

################################################################################

main()