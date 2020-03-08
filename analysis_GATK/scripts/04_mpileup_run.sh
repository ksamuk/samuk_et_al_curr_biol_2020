#!/bin/bash


#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH -J mpile_
#SBATCH -o tmp/mpile-%j.out

homedir="/dscrhome/kms173/noor2/kms173/gt_seq2"
ref="/dscrhome/kms173/noor2/kms173/iso_seq/data/ref_genome/dpse-all-chromosome-r3.04.fasta"

#bwa="~/bin/bwa"


ls -1 data/bam/*.bam | sed 's/data\/bam\///g' | sed 's/_sorted\.bam//g' > tmp/bam_list.txt

while read prefix
do

sbatch scripts/mpileup.sh $prefix

done < tmp/bam_list.txt

