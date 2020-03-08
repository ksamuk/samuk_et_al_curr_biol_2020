#!/bin/bash

#SBATCH --mem=8G
#SBATCH -p scavenger
#SBATCH --cpus-per-task=1
#SBATCH -J hap-call
#SBATCH -o tmp/hapcall/hap-call-%j.out

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
project="gt_seq_all"

# tools
bin_folder="/dscrhome/kms173/bin"
java_1_8="$bin_folder/jdk1.8.0_101/bin/java"
tabix="$bin_folder/tabix-0.2.6"
vcftools="$bin_folder/vcftools_0.1.11/bin"
picardtools="$bin_folder/picard.jar"
gatk="$bin_folder/GenomeAnalysisTK.jar"
bwa="$bin_folder/bwa"

export _JAVA_OPTIONS="-Xmx8g"
export _JAVA_OPTIONS="-XX:MaxDirectMemorySize=4G"

# the reference genome (indexed)
ref="/dscrhome/kms173/noor2/kms173/iso_seq/data/ref_genome/dpse-all-chromosome-r3.04.fasta"

mkdir -p /work/kms173/tmp
export TMPDIR=/work/kms173/tmp
export _JAVA_OPTIONS="-Djava.io.tmpdir=/work/kms173/tmp"

# convert ubams to FASTQ, pipe to BWA for alignment and then merged back raw data with aligned file
# (GATK best practices)
# markdup skipped for GT-seq


prefix=$1

	#Call GATK HaplotypeCaller
$java_1_8 -Xmx8g -Djava.io.tmpdir=/work/kms173/tmp -jar $gatk \
	-nct 1 \
	-R $ref \
	-log $log/"$prefix"_HaplotypeCaller.log \
	-T HaplotypeCaller \
	-I $bam/"$prefix".sortrg.bam \
	-L tmp/pseudo_intervals.list \
	--emitRefConfidence GVCF \
	--max_alternate_alleles 6 \
	-o $gvcf/$prefix.g.vcf
	