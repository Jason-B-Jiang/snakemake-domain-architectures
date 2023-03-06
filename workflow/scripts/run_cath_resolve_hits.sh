#!/bin/bash

# -----------------------------------------------------------------------------
#
# Get orthogroup domain architectures with cath-resolve-hits
#
# Jason Jiang - Created: Feb/17/2023
#               Last edited: Mar/03/2023

# Command line arguments taken:
#      $1 = path to hmmscan output for an orthogroup
#      $2 = path to save resolved hmmscan hits for orthogroup to
#
# Reinke Lab - Microsporidia orthologs
#
# -----------------------------------------------------------------------------

# parse in command line arguments
hmmscan_hits=$1
out=$2
orthogroup=$(basename $hmmscan_hits)

function main {
    # initialize temporary folder for storing split hmmscan hits
    mkdir -p $(pwd)/resources/tmp

    # initialize folder for storing resolved domain architectures from hmmscan
    # hits by cath-resolve-hits
    mkdir -p $out

    # run python file to create split hmmscan files in temporary folder
    python workflow/scripts/split_hmmscan_file.py ${hmmscan_hits}

    # run cath-resolve-hits on each split hmmscan file
    for split in $(pwd)/resources/tmp/${orthogroup}/*
    do
        cath-resolve-hits --input-format hmmscan_out \
            --worst-permissible-bitscore 0.1 \
            ${split} > ${out}/$(basename ${split})
    done
}

main