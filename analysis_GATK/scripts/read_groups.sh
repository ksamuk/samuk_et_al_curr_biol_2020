#!/bin/bash

#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH -J read_group
#SBATCH -o tmp/haplo-%j.out

# folders
proj_home="/dscrhome/kms173/noor2/kms173/gt_seq_rev"
work="/dscrhome/kms173/work/gt_seq_rev"
bam="$proj_home/data/bam"
log="$work/log"

# plate, barcode and project ids
project="gt_seq_run1"

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

prefix=$1

#Set MAPQ to 0 if it is unmapped, and trims alignments that overhang the reference sequence
~/bin/jdk1.8.0_101/bin/java -Djava.io.tmpdir=./tmp -jar ~/bin/picard.jar CleanSam \
 VALIDATION_STRINGENCY=LENIENT INPUT=$bam/${prefix}_sorted.bam OUTPUT=$work/data/bam/$prefix.clean.bam 2> $log/$prefix.cleansam.log
#Sort the file by the position
~/bin/jdk1.8.0_101/bin/java -Djava.io.tmpdir=./tmp -jar ~/bin/picard.jar SortSam \
VALIDATION_STRINGENCY=LENIENT INPUT=$work/data/bam/$prefix.clean.bam OUTPUT=$work/data/bam/$prefix.sort.bam SORT_ORDER=coordinate 2> $log/$prefix.sortsam.log

#Add read groups to data
~/bin/jdk1.8.0_101/bin/java -Djava.io.tmpdir=./tmp -jar ~/bin/picard.jar AddOrReplaceReadGroups \
I=$work/data/bam/$prefix.sort.bam O=$work/data/bam/$prefix.sortrg.bam SORT_ORDER=coordinate RGID=$prefix RGLB=$project RGPL=ILLUMINA RGPU=$project RGSM=$prefix CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT 2> $log/$prefix.addRG.log

	

	