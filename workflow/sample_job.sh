#!/bin/bash
#SBATCH --time=18:00:00
#SBATCH --account=*
#SBATCH --mem-per-cpu 2G
#SBATCH --cpus-per-task=5
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=*

# modules to load in on compute canada's graham server
module load python/3.10
module load singularity/3.8

source ../venv_snakemake/bin/activate

snakemake --cores all --use-singularity

deactivate
