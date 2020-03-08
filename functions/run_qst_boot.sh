#!/bin/bash
#SBATCH --mem=4G
#SBATCH -p noor
#SBATCH --cpus-per-task=1
#SBATCH -J qst_boot

Rscript qst_boot.R
