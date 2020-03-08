#!/bin/bash

#SBATCH --mem=4G
#SBATCH --cpus-per-task=2
#SBATCH -p noor
#SBATCH -J fasta_ref
#SBATCH -o tmp/build-consensus-from-list-%j.out

# 1. builds a unified vcf containing snps and indels
# 2. takes two lists, one of samples and one of gene ids with coords
# and creates consensus fasta calls in $outfolder

# single column of sample ids that match vcf sample names
samplelist=$1

# columns: chrom, pos1, pos2, fasta header string
genelist=$2

# the outputfolder with prefix
# e.g. "data/fasta/recomb_genes"
outfolder=$3

# the consensus vcf path
# eg /dscrhome/kms173/noor2/kms173/iso_seq/dpse_isoline_hardfilter_snps_indels_FILTERED.vcf.gz
consensusvcf=$4

###############################################################################
# build the unified consensus vcf
###############################################################################

# folders for building consensus vcf
proj_home="/dscrhome/kms173/noor2/kms173/iso_seq"
ref="$proj_home/data/ref_genome/dpse-all-chromosome-r3.04.fasta"
raw="$proj_home/data/raw_unzip"
demulti="$proj_home/data/demulti"
sam="$proj_home/data/sam"
bam="$proj_home/data/bam"
log="$proj_home/log"
gvcf="$proj_home/data/gvcf2"
vcf="$proj_home/data/vcf"

# plate, barcode and project ids
project="iso_seq"

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

project="iso_seq"


# re-apply hard filters, retaining indels

# build the VCF for use in generating consensus fastas

# SNPs (basic)
# $java_1_8 -Xmx32g -jar $gatk \
#     -T SelectVariants \
#     -R $ref \
#     -V $vcf/dpse_iso_seq_raw.vcf \
#     -selectType SNP \
#     -o $vcf/new/raw_snps.vcf 
# 
# Filter SNPS
# $java_1_8 -Xmx32g -jar $gatk \
#     -T VariantFiltration \
#     --disable_auto_index_creation_and_locking_when_reading_rods \
#     -R $ref \
#     -V $vcf/new/raw_snps.vcf \
#     --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
#     --filterName "BP_snp_filter" \
#     --setFilteredGtToNocall \
#     -o $vcf/new/dpse_iso_seq_hard_filtered_snps.vcf
#     
# INDELs (basic)
# $java_1_8 -Xmx32g -jar $gatk \
#     -T SelectVariants \
#     -R $ref \
#     -V $vcf/dpse_iso_seq_raw.vcf \
#     -selectType INDEL \
#     -o $vcf/new/raw_indels.vcf 
# 
# Filter INDELs
# $java_1_8 -Xmx32g -jar $gatk \
#     -T VariantFiltration \
#     --disable_auto_index_creation_and_locking_when_reading_rods \
#     -R $ref \
#     -V $vcf/new/raw_indels.vcf \
#     --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
#     --filterName "BP_indel_filter" \
#     --setFilteredGtToNocall \
#     -o $vcf/new/dpse_iso_seq_hard_filtered_indels.vcf
#     
#     
# Combine hard filtered SNPs and INDELs
# $java_1_8 -Xmx32g -jar $gatk \
#    -T CombineVariants \
#    -R $ref \
#    --variant:SNPs $vcf/new/dpse_iso_seq_hard_filtered_snps.vcf \
#    --variant:INDELs $vcf/new/dpse_iso_seq_hard_filtered_indels.vcf \
#    -o $vcf/new/dpse_isoline_hardfilter_snps_indels.vcf \
#    -genotypeMergeOptions PRIORITIZE \
#    -priority INDELs,SNPs


# compress with bgzip
# bgzip -c $vcf/dpse_isoline_hardfilter_snps_indels.vcf > dpse_iso_seq_hf_bg.vcf.gz
# tabix -p vcf dpse_iso_seq_hf_bg.vcf.gz


###############################################################################
# perform the actual consensus generation via bcftools
###############################################################################


# create consensus files

rm data/fasta_genes/*

while read geneid 
do

  # parse the gene info file for fields
  
	chrom=$(echo $geneid | awk '{print $1}')
	pos1=$(echo $geneid | awk '{print $2}')
	pos2=$(echo $geneid | awk '{print $3}')
	header=$(echo $geneid | awk '{print $4}')
	genename=$(echo $geneid | awk '{print $5}')
  cds=$(echo $geneid | awk '{print $6}')
  strand=$(echo $geneid | awk '{print $7}')
	
	while read sample
	do
	
		# get the actual sample name
		
		samplename=$(grep -w $sample meta/fasta_samplelist_recode.txt | awk '{print $2}')
	  
    # build the consensus sequence for the sample
    
		samtools faidx $ref ${chrom}:${pos1}-${pos2} | bcftools consensus $consensusvcf --iupac-codes -s $sample -o ${outfolder}_${genename}_${samplename}_single.fasta
		
    # inject the header

		echo $header | sed "s/>/>${samplename}_/g" > tmp/fasta_tmp.tmp
		sed '1d' ${outfolder}_${genename}_${samplename}_single.fasta >> tmp/fasta_tmp.tmp
		cp tmp/fasta_tmp.tmp ${outfolder}_${genename}_${samplename}_single.fasta
		rm tmp/fasta_tmp.tmp
	
	done < $samplelist
	
  # concatenate the gene files for each sample into a single file for each gene
  
	ls ${outfolder}_${genename}_*_single.fasta > tmp/gene_file_list.txt
	
	while read genefile
	do
		
		cat $genefile >> ${outfolder}_${genename}_all.fasta
		rm $genefile
	
	done < tmp/gene_file_list.txt
	
	rm tmp/gene_file_list.txt
 
  # add the cds to end of the file
  
  echo $header | sed 's/>/>CDS_/g' > tmp/fasta_tmp.tmp
  echo $cds >> tmp/fasta_tmp.tmp
  
  cat tmp/fasta_tmp.tmp >> ${outfolder}_${genename}_all.fasta
 
	
done < $genelist

