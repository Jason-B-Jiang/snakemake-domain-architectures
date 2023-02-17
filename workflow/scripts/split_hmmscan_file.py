# ------------------------------------------------------------------------------
#
# Split hmmscan file into separate files for each sequence
#
# Jason Jiang - Created: Feb/17/2023
#               Last edited: Feb/17/2023
#
# Helper script for run_cath_resolve_hits.sh
#
# Command line arguments:
#       $1 = path to hmmscan output for an orthogroup
#
# Reinke Lab - Microsporidia orthologs
#
# ------------------------------------------------------------------------------

import sys
import os
import re
from typing import List

################################################################################

def main():
    hmmscan_hits = sys.argv[1]

    orthogroup = re.search('OG\d{7}', hmmscan_hits)[0]
    os.mkdir(f"{os.getcwd()}/resources/tmp/{orthogroup}")

    with open(sys.argv[1]) as hmmscan_file:
        hmmscan_splits = hmmscan_file.read().split('\n//')
    
    for split in hmmscan_splits:
        if 'Scores for complete sequence (score includes all domains)' in split \
            and len(split) > 100:

            query_name = re.sub('(Query: +| +\[)', '', re.search('Query: +.+ +\[', split)[0])
            
            with open(f"{os.getcwd()}/resources/tmp/{orthogroup}/{query_name}", 'w') as f:
                f.write(split)

################################################################################

if __name__ == '__main__':
    main()
