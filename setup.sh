#!/bin/bash

# create virtualenv to use snakemake on cluster, and initialize other folders
# needed to run snakemake workflow
if [ ! -d "./venv_snakemake" ];
then
    # uncomment and load in the appropriate python module in your cluster
    # module load python/3.10

    virtualenv venv_snakemake
    source venv_snakemake/bin/activate

    pip install --upgrade pip
    pip install snakemake==7.22.0

    # initialize necessary folders for pipeline
    mkdir -p resources results

    # deactivate the virtual environment we just activated
    deactivate
else
    echo "run 'source venv_snakemake/bin/activate' to use Snakemake"
fi