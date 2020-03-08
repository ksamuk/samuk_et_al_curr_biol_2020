#!/bin/bash
#SBATCH -J concat
#SBATCH -p common
#SBATCH -o tmp/concat/concat-%j.out

# take a sample name and an output (root) directory
sample=$1
outdir=$2

echo "concatenating $sample..."
  
slug=$(echo ${sample} | sed 's/[0-9]*\///g')

echo "with slug $slug..."
  
cat $outdir/data/demulti/lane1/$sample \
  $outdir/data/demulti/lane2/$sample \
  $outdir/data/demulti/lane3/$sample \
  $outdir/data/demulti/lane4/$sample > /hpchome/noor/kms173/work/gt_seq_rev/data/demulti/$slug
  

