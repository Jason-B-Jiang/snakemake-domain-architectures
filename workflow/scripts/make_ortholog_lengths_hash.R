# -----------------------------------------------------------------------------
#
# Create hashmap of ortholog protein lengths
#
# Jason Jiang - Created: Mar/09/2023
#               Last edited: Mar/09/2023
#
# Reinke Lab - Microsporidia orthologs
#
# -----------------------------------------------------------------------------

suppressMessages(library(tidyverse))
suppressMessages(library(seqinr))

################################################################################

main <- function() {
  # ---------------------------------------------------------------------------
  # Command line arguments:
  #   $1 = filepath to folder with species proteomes (fasta files)
  #   $2 = filepath to save ortholog length hashmap to
  # ---------------------------------------------------------------------------
  args <- commandArgs(trailingOnly = TRUE)
  proteomes <- args[1]
  out <- args[2]
  
  write_rds(make_orthologs_length_hash(proteomes), out)
}

################################################################################

## Helper functions

make_orthologs_length_hash <- function(proteomes) {
  # ---------------------------------------------------------------------------
  # Docstring goes here
  # ---------------------------------------------------------------------------
  species_lengths <- new.env()
  
  for (proteome in list.files(proteomes, full.names = TRUE)) {
    ortholog_lengths <- new.env()
    species <- str_remove(basename(proteome), '\\.fa')
    
    for (ortholog in read.fasta(proteome, seqtype = 'AA')) {
      ortholog_lengths[[getName(ortholog)]] <- getLength(ortholog)
    }
    
    species_lengths[[species]] <- ortholog_lengths
  }
  
  return(species_lengths)
}

################################################################################

main()