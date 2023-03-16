# -----------------------------------------------------------------------------
#
# Get alignment ends for local alignments between single-copy ortholog pairs
#
# Jason Jiang - Created: 2023/03/16
#               Last edited: 2023/03/16
#
# Reinke Lab - Microsporidia Orthologs Project
#
#
# -----------------------------------------------------------------------------

suppressMessages(library(tidyverse))
suppressMessages(library(seqinr))
suppressMessages(library(Biostrings))

################################################################################

# Parameters for Smith-Waterman alignment
GAP_OPEN = 11.0
GAP_EXTEND = 1.0
ALIGNMENT_TYPE = 'local'

################################################################################

main <- function() {
  # ---------------------------------------------------------------------------
  # Command line arguments:
  #   $1 = folder w/ fasta sequences for all orthogroups
  #   $2 = all single-copy ortholog pairs with yeast from orthofinder
  #   $3 = output filepath
  # ---------------------------------------------------------------------------
  args <- commandArgs(trailingOnly = T)
  
  orthogroup_seqs <- args[1]
  
  orthogroups <- read_csv(args[2], show_col_types = F) %>%
    select(orthogroup, species, is_microsp, yeast_ortholog, species_ortholog,
           yeast_len, species_len)
  
  out <- args[3]
  
  # make hashtable mapping all orthologs in orthogroups to their sequences
  orthogroup_seqs_hash <- make_orthogroup_seqs_hash(orthogroup_seqs,
                                                    unique(orthogroups$orthogroup))
  
  # align species-yeast ortholog pairs with local alignment and get the end
  # of the species alignment
  orthogroups <- orthogroups %>%
    rowwise() %>%
    mutate(alignment_end = get_c_term_end(
      orthogroup_seqs_hash[[orthogroup]][[species_ortholog]],
      orthogroup_seqs_hash[[orthogroup]][[yeast_ortholog]],
      is_microsp,
      species_len,
      yeast_len
    ))
  
  write_csv(orthogroups, out)
}
  
################################################################################

## Helper functions

make_orthogroup_seqs_hash <- function(orthogroup_seqs,
                                      orthogroups_of_interest) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  base_dir <- dirname(list.files(orthogroup_seqs, full.names = T)[1])
  seq_fasta <- unname(sapply(orthogroups_of_interest,
                             function(og) {str_c(base_dir, '/', og, '.fa')}))
  
  # TODO - vectorize this loop
  orthogroup_seqs_hash <- new.env()
  for (seq in seq_fasta) {
    orthogroup_seqs <- new.env()
    fa <- read.fasta(seq, seqtype = 'AA', as.string = T)
    orthologs <- getName(fa)
    
    for (ortholog in orthologs) {
      orthogroup_seqs[[ortholog]] <- fa[[ortholog]][1]
    }
    
    orthogroup_seqs_hash[[basename(str_remove(seq, '\\.fa'))]] <-
      orthogroup_seqs
  }
  
  return(orthogroup_seqs_hash)
}


get_c_term_end <- function(species_seq, yeast_seq, is_microsp, species_len,
                           yeast_len) {
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  if (species_len >= yeast_len) {
    return(NA)
  }
  
  if (is_microsp) {
    MATRIX = 'BLOSUM45'
  } else {
    MATRIX = 'BLOSUM62'
  }
  
  local_alignment <- pairwiseAlignment(species_seq, yeast_seq,
                                       type = ALIGNMENT_TYPE,
                                       substitutionMatrix = MATRIX,
                                       gapOpening = GAP_OPEN,
                                       gapExtension = GAP_EXTEND)
  
  # return residue in species ortholog where alignment ends, indicating beginning
  # of c-terminus truncation
  return(end(subject(local_alignment)))
}

################################################################################

main()

################################################################################

ogs <- read_csv('single_copy_orthogroups.csv')
ref_species <- 'E_brev'
outgroup <- 'E_hell'
ortho_lengths <- read_rds('../resources/ortholog_lengths.rds')

species_orthos <- select(ogs, -ref_species) %>% rename(Orthogroup = 'orthogroup')
ref_orthos <- select(ogs, Orthogroup, ref_species) %>% rename(Orthogroup = 'orthogroup')
colnames(ref_orthos)[2] = 'ref_ortholog'

species_orthos <- species_orthos %>%
  pivot_longer(cols = colnames(species_orthos)[2 : ncol(species_orthos)],
               names_to = "species",
               values_to = "ortholog") %>%
  full_join(ref_orthos, by = 'orthogroup') %>%
  # remove entries with missing ortholog for species, or non-single copy ortholog
  filter(!is.na(ortholog), !str_detect(ortholog, ',')) %>%
  mutate(is_outgroup = species %in% outgroup) %>%
  rowwise() %>%
  mutate(ortholog_length = ortho_lengths[[orthogroup]][[ortholog]],
         ref_ortholog_length = ortho_lengths[[orthogroup]][[ref_ortholog]])