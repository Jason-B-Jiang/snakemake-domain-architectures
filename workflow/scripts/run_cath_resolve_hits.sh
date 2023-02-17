#!/bin/bash

# -----------------------------------------------------------------------------
#
# Get orthogroup domain architectures with cath-resolve-hits
#
# Jason Jiang - Created: Feb/17/2023
#               Last edited: Feb/17/2023

# Command line arguments taken:
#      $1 = path to hmmscan output for an orthogroup
#      $2 = path to save resolved hmmscan hits for orthogroup to
#
# Reinke Lab - Microsporidia orthologs
#
# -----------------------------------------------------------------------------

# check working directory of script at runtime
# load in command line arguments (input, output)
# run python helper script to split hmmscan files into resources/tmp folder
# for each split hmmscan file in resources/tmp/orthogroup, run crh
# remove resources/tmp folder
# done

# parse in command line arguments
hmmscan_hits=$1
out=$2

function main {
    # initialize temporary folder for storing split hmmscan hits
    mkdir -p $(pwd)/resources/tmp

    # run python file to create split hmmscan files in temporary folder
    python workflow/scripts/split_hmmscan_file.py ${hmmscan_hits}

    # rm -r "$PWD"/resources/tmp
}

main