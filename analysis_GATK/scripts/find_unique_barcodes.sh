#!/bin/bash


#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH -J demulti_run
#SBATCH -o tmp/demulti-%j.out

homedir="/dscrhome/kms173/noor2/kms173/gt_seq_rev"
outdir="/dscrhome/kms173/work/gt_seq_rev"
raw="/dscrhome/kms173/noor2/kms173/gt_seq_rev/data/raw/Samuk_5087_180803A5/"



#grep -A1 @ $raw/5087_pool_S1_L004_I2_001.fastq | grep -v ^- | grep -v @ | sort | uniq > barcodes_test.txt
#grep -A1 @ $raw/5087_pool_S1_L004_I2_001.fastq | grep -v ^- | grep -v @ | sort > ~/work/all_barcodes_test.txt
uniq -c ~/work/all_barcodes_test.txt > barcodes_counts.txt
