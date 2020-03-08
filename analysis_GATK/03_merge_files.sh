#!/bin/bash


#SBATCH --mem=4G
#SBATCH --cpus-per-task=2
#SBATCH -J concat_
#SBATCH -o concat-%j.out

homedir="/dscrhome/kms173/noor2/kms173/gt_seq_rev"
outdir="/hpchome/noor/kms173/work/gt_seq_rev"
raw="/hpchome/noor/kms173/work/gt_seq_rev/run2/data/raw"

# RUN1


find $outdir/run1/data/demulti -name "*.fq" | sed 's/\/hpchome\/noor\/kms173\/work\/gt_seq_rev\/run1\/data\/demulti\/lane[0-9]\///g' | sort | uniq > tmp/run1_samples.txt
find $outdir/run2/data/demulti -name "*.fq" | sed 's/\/hpchome\/noor\/kms173\/work\/gt_seq_rev\/run2\/data\/demulti\/lane[0-9]\///g' | sort | uniq > tmp/run2_samples.txt

while read sample
do

sbatch scripts/concatenate.sh $sample $outdir/run1

done < tmp/run1_samples.txt



while read sample
do

sbatch scripts/concatenate.sh $sample $outdir/run2

done < tmp/run2_samples.txt