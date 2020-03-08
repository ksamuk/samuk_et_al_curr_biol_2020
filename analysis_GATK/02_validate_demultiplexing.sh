#!/bin/bash


#SBATCH --mem=6G
#SBATCH -p noor
#SBATCH --cpus-per-task=1
#SBATCH -J demulti_run
#SBATCH -o tmp/demulti-launch-%j.out

homedir="/dscrhome/kms173/noor2/kms173/gt_seq_rev"
outdir1="/hpchome/noor/kms173/work/gt_seq_rev/run1"
outdir2="/hpchome/noor/kms173/work/gt_seq_rev/run2"


#find /hpchome/noor/kms173/work/gt_seq_rev/run1 -name '*.fq' -exec ls -lh {} \; > tmp/demulti_files_run1.list
#find /hpchome/noor/kms173/work/gt_seq_rev/run2 -name '*.fq' -exec ls -lh {} \; > tmp/demulti_files_run2.list

find $outdir1 -name '*.fq' -exec ls -lh {} \; > tmp/demulti_files_run1.list
find $outdir2 -name '*.fq' -exec ls -lh {} \; > tmp/demulti_files_run2.list

