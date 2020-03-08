#!/bin/bash

#SBATCH --mem=8G
#SBATCH -p common
#SBATCH --cpus-per-task=1
#SBATCH -J demulti
#SBATCH -o tmp/demulti2/demulti-%j.out

homedir="/hpc/group/noor/kms173/noor2/gt_seq_rev"
#outdir="/hpchome/noor/kms173/work/gt_seq_rev/run2"
raw="/hpc/home/kms173/work/gt_seq_rev/run${4}/data/raw"

# 3rd argument is the output(root) directory
outdir=$3
run=$4
slug=$5

conda activate genomics

# function for demultiplexing ONE plate
demult(){
  echo "run${run}, lane ${2}, i7 ${1}";
  mkdir $outdir/data/demulti/lane${2};
  mkdir $outdir/data/demulti/lane${2}/$1;
  echo "${homedir}/meta/run${run}/multx_indexes/run${run}_multx_indexes_"$1.txt;
  
  fastq-multx -d 2 -m 0 -B "${homedir}/meta/run${run}/multx_indexes/run${run}_multx_indexes_"$1.txt $raw/${slug}_S1_L00${2}_I1_001.fastq $raw/${slug}_S1_L00${2}_I2_001.fastq $raw/${slug}_S1_L00${2}_R1_001.fastq $raw/${slug}_S1_L00${2}_R2_001.fastq -o n/a -o n/a -o $outdir/data/demulti/lane$2/$1/r1.%.fq -o $outdir/data/demulti/lane$2/$1/r2.%.fq;
  
  # have to remove these for per plate, otherwise we end up with 51x the raw data...
  rm $outdir/data/demulti/lane$2/$1/r2.unmatched.fq;
  rm $outdir/data/demulti/lane$2/$1/r1.unmatched.fq;
}


demult $1 $2 $outdir $run $slug

echo "COMPLETED SUCCESSFULLY"



