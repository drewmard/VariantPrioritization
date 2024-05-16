#!/bin/bash

# # sbatch --account=smontgom --partition=batch --time=1-1:00:00 --mem=64G --nodes=1 --ntasks=1 /oak/stanford/groups/smontgom/amarder/neuro-variants/scripts/snps/finemap/run_finemap_pipeline.sh Alzheimers_Bellenguez_2022 1

# module load R/4.1.2
mamba init
mamba activate r

trait=$1
startInd=$2

path_to_script=/oak/stanford/groups/smontgom/amarder/neuro-variants/scripts/snps/finemap/finemap_pipeline.R
Rscript $path_to_script $trait $startInd
