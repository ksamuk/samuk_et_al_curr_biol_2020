#!/bin/bash

#SBATCH --mem=220000
#SBATCH --cpus-per-task=13
#SBATCH -p noor
#SBATCH -J il_hmlgy
#SBATCH -o tmp/il_hmlgy-%j.out

source activate r_env

Rscript scripts/07_estimate_homology_illumina.R
