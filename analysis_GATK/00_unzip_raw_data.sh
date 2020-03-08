#!/bin/bash

#SBATCH -p noor
#SBATCH -J unzip

mkdir ~/work/gt_seq_rev
mkdir ~/work/gt_seq_rev/run1
mkdir ~/work/gt_seq_rev/run1/data
mkdir ~/work/gt_seq_rev/run1/data/raw

mkdir ~/work/gt_seq_rev/run2
mkdir ~/work/gt_seq_rev/run2/data
mkdir ~/work/gt_seq_rev/run2/data/raw


#cd data/raw/Samuk_5228_181025A5

#for a in *.gz; do gunzip -c $a > ~/work/gt_seq_rev/run2/data/raw/`echo $a | sed s/.gz//`; done


cd /datacommons/noor2/kms173/gt_seq_rev/data/raw/Samuk_5087_180803A5

for a in *fastq.gz; do gunzip -c $a > ~/work/gt_seq_rev/run1/data/raw/`echo $a | sed s/.gz//`; done


