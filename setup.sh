#!/bin/bash

# create virtualenv to use snakemake on cluster, and initialize other folders
# needed to run snakemake workflow
if [ ! -d "./snakemake_env" ];
then
    # uncomment and load in the appropriate python module in your cluster
    # module load python/3.10

    virtualenv snakemake_env
    source snakemake_env/bin/activate

    pip install --upgrade pip
    pip install snakemake==7.22.0

    # initialize necessary folders for pipeline
    mkdir -p resources results
else
    echo "run 'source snakemake_env/bin/activate' to use Snakemake"
fi