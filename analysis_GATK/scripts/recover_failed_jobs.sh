#!/bin/bash

#SBATCH --mem=4G
#SBATCH --cpus-per-task=2
#SBATCH -J align_fix
#SBATCH -p scavenger
#SBATCH -o tmp/align-%j.out

homedir="/dscrhome/kms173/work/gt_seq_rev"
ref="/dscrhome/kms173/noor2/kms173/iso_seq/data/ref_genome/dpse-all-chromosome-r3.04.fasta"

sacct -u kms173 --format=User,JobID,Jobname,Partition,state,start,end,elapsed | grep scavenger | grep FAILED | grep align | grep -v 2018-08-10T11 | awk '{ print $2 }' > failed_jobs.txt


while read job_id; do

prefix=$(cat tmp/align-$job_id.out | head -n1 | sed 's/aligning //g' | sed 's/\.\.\.//g')

sbatch scripts/align.sh $prefix

done < failed_jobs.txt
