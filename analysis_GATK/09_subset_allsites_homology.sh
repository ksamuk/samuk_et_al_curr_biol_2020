#!/bin/bash

#SBATCH --cpus-per-task=2
#SBATCH -p noor
#SBATCH -J regions
#SBATCH -o tmp/regions-%j.out

mkdir data/allsites_recomb

while read interval
do

echo $interval > tmp/interval_tmp.txt

siteshort=$(echo $interval | sed 's/ /_/g')

echo $siteshort

tabix -h -R tmp/interval_tmp.txt ~/work/iso_seq/iso_seq_all_sites.vcf.gz | gzip > data/allsites_recomb/${siteshort}_recomb_regions.vcf.gz

done < meta/recom_intervals.txt

#tabix -h -R meta/recom_intervals.txt ~/work/iso_seq/iso_seq_all_sites.vcf.gz | gzip > data/iso_seq_all_sites_recomb_regions.vcf.gz
