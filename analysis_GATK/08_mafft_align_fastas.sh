#!/bin/bash

#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH -J mafft
#SBATCH -p noor
#SBATCH -o tmp/fasta_align-%j.out

# folders
proj_home="/dscrhome/kms173/work/gt_seq_rev"
ref="$proj_home/data/ref_genome/dpse-all-chromosome-r3.04.fasta"
raw="$proj_home/data/raw_unzip"
demulti="$proj_home/data/demulti"
sam="$proj_home/data/sam"
bam="$proj_home/data/bam"
log="$proj_home/log"
gvcf="$proj_home/data/gvcf"

# plate, barcode and project ids
project="gt_seq_test2"

# tools
bin_folder="/dscrhome/kms173/bin"
java_1_8="$bin_folder/jdk1.8.0_101/bin/java"
tabix="$bin_folder/tabix-0.2.6"
vcftools="$bin_folder/vcftools_0.1.11/bin"
picardtools="$bin_folder/picard.jar"
gatk="$bin_folder/GenomeAnalysisTK.jar"
bwa="$bin_folder/bwa"


ls data/fasta_genes/*.fasta | sed "s/.*\///g" | sed 's/.fasta//g' > tmp/fasta_gene_list.txt

mkdir data/fasta_genes/aligned

rm data/fasta_genes/aligned/*

while read fastafile
do

/dscrhome/kms173/bin/mafft-linux64/mafft.bat --thread 4 --maxiterate 1000 data/fasta_genes/$fastafile.fasta > data/fasta_genes/aligned/${fastafile}_aligned.fasta

done < tmp/fasta_gene_list.txt

