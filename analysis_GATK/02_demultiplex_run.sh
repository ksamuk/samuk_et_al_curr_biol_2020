#!/bin/bash


#SBATCH --mem=6G
#SBATCH -p noor
#SBATCH --cpus-per-task=1
#SBATCH -J demulti_run
#SBATCH -o tmp/demulti-launch-%j.out

homedir="/hpc/group/noor/kms173/noor2/gt_seq_rev"
outdir1="/hpc/home/kms173/work/gt_seq_rev/run1"
outdir2="/hpc/home/kms173/work/gt_seq_rev/run2"

echo "preparing to launch jobs..."

# LANE 1

mkdir $outdir1/data/demulti

rm -r $outdir1/data/demulti/lane1
mkdir $outdir1/data/demulti/lane1

rm -r $outdir1/data/demulti/lane2
mkdir $outdir1/data/demulti/lane2

rm -r $outdir1/data/demulti/lane3
mkdir $outdir1/data/demulti/lane3

rm -r $outdir1/data/demulti/lane4
mkdir $outdir1/data/demulti/lane4

# test

#sbatch scripts/01_demultiplex_one_lane.sh 1 1 /hpc/home/noor/kms173/work/gt_seq_rev/run1 1 5087_pool

# lane 1
for i in `seq 1 49`; do
  sbatch scripts/01_demultiplex_one_lane.sh $i 1 $outdir1 1 5087_pool
done

# lane 1
for i in `seq 1 49`; do
  sbatch scripts/01_demultiplex_one_lane.sh $i 2 $outdir1 1 5087_pool
done

# lane 1
for i in `seq 1 49`; do
  sbatch scripts/01_demultiplex_one_lane.sh $i 3 $outdir1 1 5087_pool
done

# lane 1
for i in `seq 1 49`; do
  sbatch scripts/01_demultiplex_one_lane.sh $i 4 $outdir1 1 5087_pool
done 


# top up i7 51
sbatch scripts/01_demultiplex_one_lane.sh 51 1 $outdir1 1 5087_pool
sbatch scripts/01_demultiplex_one_lane.sh 51 2 $outdir1 1 5087_pool
sbatch scripts/01_demultiplex_one_lane.sh 51 3 $outdir1 1 5087_pool
sbatch scripts/01_demultiplex_one_lane.sh 51 4 $outdir1 1 5087_pool


# LANE 2

mkdir $outdir2/data/demulti

rm -r $outdir2/data/demulti/lane1
mkdir $outdir2/data/demulti/lane1

rm -r $outdir2/data/demulti/lane2
mkdir $outdir2/data/demulti/lane2

rm -r $outdir2/data/demulti/lane3
mkdir $outdir2/data/demulti/lane3

rm -r $outdir2/data/demulti/lane4
mkdir $outdir2/data/demulti/lane4


for i in `seq 1 51`; do
  sbatch scripts/01_demultiplex_one_lane.sh $i 1 $outdir2 2 5228_pool
done


for i in `seq 1 51`; do
  sbatch scripts/01_demultiplex_one_lane.sh $i 2 $outdir2 2 5228_pool
done


for i in `seq 1 51`; do
  sbatch scripts/01_demultiplex_one_lane.sh $i 3 $outdir2 2 5228_pool
done


for i in `seq 1 51`; do
  sbatch scripts/01_demultiplex_one_lane.sh $i 4 $outdir2 2 5228_pool
done 

