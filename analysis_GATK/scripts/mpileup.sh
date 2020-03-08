#!/bin/bash


#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH -J mpile
#SBATCH -o tmp/mpile-%j.out

homedir="/dscrhome/kms173/noor2/kms173/gt_seq2"
ref="/dscrhome/kms173/noor2/kms173/iso_seq/data/ref_genome/dpse-all-chromosome-r3.04.fasta"

bamfile=$1

samtools mpileup -d 100 -R -uf $ref data/bam/${bamfile}_sorted.bam > data/bcf/${bamfile}.bcf  
~/bin/bcftools call --threads 2 --multiallelic-caller data/bcf/${bamfile}.bcf   -o data/vcf/${bamfile}.vcf 
rm data/bcf/${bamfile}.bcf
cd data/vcf
bgzip ${bamfile}.vcf 
tabix -p vcf ${bamfile}.vcf.gz