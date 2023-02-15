#!/bin/bash

# create virtualenv to use snakemake on cluster, and initialize other folders
# needed to run snakemake workflow
if [ ! -d "./snakemake_env" ];
then
    module load python/3.10

    virtualenv --no-download snakemake_env
    source snakemake_env/bin/activate

    pip install --no-index --upgrade pip
    pip install snakemake

    # initialize necessary folders for pipeline
    mkdir -p resources results
else
    echo "run 'source snakemake_env/bin/activate' to use Snakemake"
fi