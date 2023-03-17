# -----------------------------------------------------------------------------
#
# Get c-terminus truncation between orthologs and reference orthologs
#
# Jason Jiang - Created: 2023/03/16
#               Last edited: 2023/03/17
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
  orthogroups <- read_csv(args[2], show_col_types = F)
  out <- args[3]
  
  orthogroup_seqs <- '../../results/OrthoFinder/Results_OrthoFinder/Orthogroup_Sequences'
  orthogroups <- read_csv('../../results/single_copy_orthogroups.csv')
  out <- '../../results/temp.csv'
  
  # make hashtable mapping all orthologs in orthogroups to their sequences
  orthogroup_seqs_hash <- make_orthogroup_seqs_hash(orthogroup_seqs,
                                                    unique(orthogroups$orthogroup))
  
  # align species-yeast ortholog pairs with local alignment and get the end
  # of the species alignment
  orthogroups <- orthogroups %>%
    rowwise() %>%
    mutate(alignment_end = get_c_term_end(
      orthogroup_seqs_hash[[orthogroup]][[ortholog]],
      orthogroup_seqs_hash[[orthogroup]][[ref_ortholog]],
      is_outgroup,
      ortholog_length,
      ref_ortholog_length
    ))
  
  write_csv(orthogroups, out)
}
  
################################################################################

## Helper functions

# Functions for doing local alignments between orthologs and reference orthologs,
# and getting extent of c-terminus truncation for the ortholog

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


get_c_term_end <- function(ortho_seq, ref_seq, is_outgroup, ortho_len,
                           ref_len) {
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  if (ortho_len >= ref_len) {
    return(NA)
  }
  
  if (!is_outgroup) {
    MATRIX = 'BLOSUM45'
  } else {
    MATRIX = 'BLOSUM62'
  }
  
  local_alignment <- pairwiseAlignment(ortho_seq, ref_seq,
                                       type = ALIGNMENT_TYPE,
                                       substitutionMatrix = MATRIX,
                                       gapOpening = GAP_OPEN,
                                       gapExtension = GAP_EXTEND)
  
  # return residue in species ortholog where alignment ends, indicating beginning
  # of c-terminus truncation in species ortholog, relative to reference ortholog
  return(end(subject(local_alignment)))
}

# Functions for classifying lost c-terminus residues in orthologs

################################################################################

main()