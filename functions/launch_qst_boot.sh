#!/bin/bash
#SBATCH --mem=4G
#SBATCH -p noor
#SBATCH --cpus-per-task=1
#SBATCH -J qst_

for i in {1..100}
do
  sbatch run_qst_boot.sh
done
