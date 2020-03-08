#!/bin/bash

#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH -J haplo_
#SBATCH -o tmp/haplo-%j.out

# folders
proj_home="/dscrhome/kms173/work/gt_seq_rev"
ref="$proj_home/data/ref_genome/dpse-all-chromosome-r3.04.fasta"
raw="$proj_home/data/raw_unzip"
demulti="$proj_home/data/demulti"
sam="$proj_home/data/sam"
bam="/dscrhome/kms173/work/gt_seq_rev/data/bam"
log="$proj_home/log"
gvcf="$proj_home/data/gvcf"

# plate, barcode and project ids
project="gt_seq_all"

# tools
bin_folder="/dscrhome/kms173/bin"
java_1_8="$bin_folder/jdk1.8.0_101/bin/java"
tabix="$bin_folder/tabix-0.2.6"
vcftools="$bin_folder/vcftools_0.1.11/bin"
picardtools="$bin_folder/picard.jar"
gatk="$bin_folder/GenomeAnalysisTK.jar"
bwa="$bin_folder/bwa"

export _JAVA_OPTIONS="-Xmx20g"
export _JAVA_OPTIONS="-XX:MaxDirectMemorySize=4G"

# the reference genome (indexed)
ref="/dscrhome/kms173/noor2/kms173/iso_seq/data/ref_genome/dpse-all-chromosome-r3.04.fasta"

ls -1 $bam/*.sortrg.bam | sed 's/.*\///g' | sed 's/\..*//g' > tmp/bam_list.txt

while read prefix; do

	sbatch scripts/haplotype_caller.sh $prefix
	
done < tmp/bam_list.txt
	