#!/bin/bash

#SBATCH --cpus-per-task=2
#SBATCH -p noor
#SBATCH -J homlgy
#SBATCH -o tmp/homlgy-%j.out

# create a temporary fasta with corrected headers

sed 's/|arrow//g' ~/work/iso_pac_bio/data/consensus/${1}_consensus.fasta > data/consensus_sub/${1}_temp.fasta

# subset the consensus file for the bed regions
bedtools getfasta -fi data/consensus_sub/${1}_temp.fasta \
-bed meta/recom_intervals.txt > data/consensus_sub/${1}_recomb_regions.fasta

rm data/consensus_sub/${1}_temp.*