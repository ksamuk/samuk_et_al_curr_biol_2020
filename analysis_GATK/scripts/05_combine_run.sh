#!/bin/bash


#SBATCH --mem=20G
#SBATCH --cpus-per-task=4
#SBATCH -J vcf-merge
#SBATCH -o tmp/vcf-merge-%j.out

homedir="/dscrhome/kms173/noor2/kms173/gt_seq2"
ref="/dscrhome/kms173/noor2/kms173/iso_seq/data/ref_genome/dpse-all-chromosome-r3.04.fasta"

ls data/vcf | sed 's/_[A-H]*[0-9]*.vcf.*//g' | sort | uniq > tmp/plate_list.txt

mkdir data/vcf_combined

while read plate; do

	sbatch scripts/combine.sh $plate
	
done < tmp/plate_list.txt