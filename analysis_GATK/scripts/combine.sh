#!/bin/bash


#SBATCH --mem=10G
#SBATCH --cpus-per-task=2
#SBATCH -J vcf-merge
#SBATCH -o tmp/vcf-merge-%j.out

homedir="/dscrhome/kms173/noor2/kms173/gt_seq2"
ref="/dscrhome/kms173/noor2/kms173/iso_seq/data/ref_genome/dpse-all-chromosome-r3.04.fasta"

plate=$1

ls data/vcf/*.vcf.gz | grep $plate | xargs vcf-merge -t -d > data/vcf_combined/gt_seq_run1_$plate.vcf