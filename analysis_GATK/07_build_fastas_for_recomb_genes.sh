#!/bin/bash


#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH -J fasta_ref
#SBATCH -p noor
#SBATCH -o tmp/fasta_ref-%j.out

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


sh scripts/create_fasta_seqs_from_list.sh \
meta/fasta_samplelist.txt \
meta/anderson_gene_locations.txt \
data/fasta_genes/recomb_genes \
/dscrhome/kms173/noor2/kms173/iso_seq/dpse_isoline_hardfilter_snps_indels_FILTERED.vcf.gz

