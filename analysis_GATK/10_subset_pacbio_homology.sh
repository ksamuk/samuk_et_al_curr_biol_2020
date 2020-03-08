#!/bin/bash

#SBATCH --cpus-per-task=2
#SBATCH -p noor
#SBATCH -J homlgy
#SBATCH -o tmp/homlgy-%j.out


ls ~/work/iso_pac_bio/data/consensus/*.fasta | sed 's/.*\///g' | sed 's/_.*//g' > consensus_list.txt

mkdir data/consensus_sub

while read consensus
do

sbatch scripts/subset_pacbio_homology.sh $consensus

done < consensus_list.txt

rm consensus_list.txt
