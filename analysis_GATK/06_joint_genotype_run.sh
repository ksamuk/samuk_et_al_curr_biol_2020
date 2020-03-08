#!/bin/bash


#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH -J genocall_
#SBATCH -p noor
#SBATCH -o tmp/genocall-%j.out

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

export _JAVA_OPTIONS="-Xmx20g"
export _JAVA_OPTIONS="-XX:MaxDirectMemorySize=4G"

# the reference genome (indexed)
ref="/dscrhome/kms173/noor2/kms173/iso_seq/data/ref_genome/dpse-all-chromosome-r3.04.fasta"

ls $proj_home/data/gvcf | sed 's/_[A-H]*[0-9]*.g.vcf.*//g' | sort | uniq | sed 's/_[0-9]*$//g' | sort | uniq > tmp/line_list.txt


while read line; do

	sbatch scripts/joint_genotype.sh $line
	
done < tmp/line_list.txt


