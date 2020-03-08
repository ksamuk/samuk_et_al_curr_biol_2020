#!/bin/bash

#SBATCH --mem=4G
#SBATCH --cpus-per-task=2
#SBATCH -J align_
#SBATCH -p noor
#SBATCH -o tmp/align-launch-%j.out

homedir="/hpchome/noor/kms173/work/gt_seq_rev/"
bam="/hpchome/noor/kms173/work/gt_seq_rev/data/bam"
demulti="/hpchome/noor/kms173/work/gt_seq_rev/data/demulti"
ref="/dscrhome/kms173/noor2/kms173/iso_seq/data/ref_genome/dpse-all-chromosome-r3.04.fasta"

mkdir $homedir/log

#bwa="~/bin/bwa"

ls $homedir/data/demulti/ | sed 's/r[12]\.//g' | sed 's/\.fq//g' | sort | uniq  > tmp/sample_list.txt

#ls ~/work/gt_seq_rev/data/demulti2 | sed 's/r[12]\.//g' | sed 's/\.fq//g' | sort | uniq  > tmp/sample_list.txt

while read prefix
do

sbatch scripts/align.sh $prefix

done < tmp/sample_list.txt
