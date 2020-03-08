#!/bin/bash

#SBATCH -p noor
#SBATCH --mem=40G
#SBATCH --cpus-per-task=4
#SBATCH -J joint-geno
#SBATCH -o tmp/joint-geno-%j.out

# folders
proj_home="/dscrhome/kms173/noor2/kms173/gt_seq_rev"
work_home="/dscrhome/kms173/work/gt_seq_rev"
ref="$proj_home/data/ref_genome/dpse-all-chromosome-r3.04.fasta"
raw="$proj_home/data/raw_unzip"
demulti="$proj_home/data/demulti"
sam="$proj_home/data/sam"
bam="$proj_home/data/bam"
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

export _JAVA_OPTIONS="-Xmx12g"
export _JAVA_OPTIONS="-XX:MaxDirectMemorySize=4G"

# the reference genome (indexed)
ref="/dscrhome/kms173/noor2/kms173/iso_seq/data/ref_genome/dpse-all-chromosome-r3.04.fasta"

export _JAVA_OPTIONS="-Djava.io.tmpdir=/work/kms173/tmp"
export TMPDIR=/work/kms173/tmp

line=$1

#ls $work_home/data/gvcf/*.g.vcf | grep $line > tmp/$line.gvcf.list

ls $work_home/data/gvcf/*.g.vcf | sort | uniq | grep $line > tmp/$line.gvcf.list

#Call GATK joint genotyper

$java_1_8 -Xmx40g -Djava.io.tmpdir=/work/kms173/tmp -jar $gatk \
  -T GenotypeGVCFs \
  -V tmp/$line.gvcf.list \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  -nt 4 \
  -L tmp/pseudo_intervals.list \
  -R $ref \
  -log $log/GenotypeNEGWGSGVCFs.log \
  -l INFO \
  -o data/gatk_vcf/$line.vcf
