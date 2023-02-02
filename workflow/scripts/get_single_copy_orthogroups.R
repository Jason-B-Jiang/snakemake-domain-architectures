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

library(tidyverse)

################################################################################

main <- function() {
  # ---------------------------------------------------------------------------
  # Command line arguments:
  # $1 = filepath to directory with OrthoFinder output
  # $2 = name of reference species to find single-copy orthogroups for
  # $3 = filepath to save single-copy orthogroups to
  # $4 = filepath to save NAMES of single-copy orthogroups to
  #      (this will help later when choosing orthogroups to run hmmscan on)
  # ---------------------------------------------------------------------------
  args <- commandArgs(trailingOnly = TRUE)
  orthogroups <- read_tsv(str_c(args[1], '/Results_OrthoFinder/Orthogroups/Orthogroups.tsv'),
                          show_col_types = FALSE)
  ref_species <- args[2]
  orthogroup_out <- args[3]
  orthogroup_names_out <- args[4]
  
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
  
  # write orthogroups as csv file to specified out filepath
  write_csv(orthogroups, orthogroup_out)
  
  # write text file of orthogroups fitting this criterion to resources folder
  # (will need when picking which orthogroups to assign domains to with hmmscan)
  write_lines(orthogroups$Orthogroup, orthogroup_names_out)
}

################################################################################

main()