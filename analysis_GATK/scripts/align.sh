#!/bin/bash

#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH -J align
#SBATCH -p common
#SBATCH -o tmp/align/align-%j.out

#homedir="/dscrhome/kms173/noor2/kms173/gt_seq_rev"

homedir="/hpchome/noor/kms173/work/gt_seq_rev"
bam="/hpchome/noor/kms173/work/gt_seq_rev/data/bam"
demulti="/hpchome/noor/kms173/work/gt_seq_rev/data/demulti"
ref="/dscrhome/kms173/noor2/kms173/iso_seq/data/ref_genome/dpse-all-chromosome-r3.04.fasta"

log="${homedir}/log"
project="gt_seq_combined"

#bwa="~/bin/bwa"
prefix=$1

echo "aligning ${prefix}..."
# align with bwa
bwa mem -t 1 $ref $homedir/data/demulti/r1.${prefix}.fq $homedir/data/demulti/r2.${prefix}.fq > $homedir/data/bam/${prefix}.sam 

#bwa mem -t 1 /dscrhome/kms173/noor2/kms173/iso_seq/data/ref_genome/dpse-all-chromosome-r3.04.fasta ~/work/gt_seq_rev/data/demulti2/r1.AFC_14_15_A01.fq ~/work/gt_seq_rev/data/demulti2/r2.AFC_14_15_A01.fq > ~/work/gt_seq_rev/data/bam/AFC_14_15_A01.sam 

# sam to bam, sort
samtools view -Sb $homedir/data/bam/${prefix}.sam  > $homedir/data/bam/${prefix}.bam
samtools sort $homedir/data/bam/${prefix}.bam > $homedir/data/bam/${prefix}_sorted.bam
samtools index $homedir/data/bam/${prefix}_sorted.bam

# remove intermediates
rm $homedir/data/bam/${prefix}.sam
rm $homedir/data/bam/${prefix}.bam

#Set MAPQ to 0 if it is unmapped, and trims alignments that overhang the reference sequence
~/bin/jdk1.8.0_101/bin/java -Djava.io.tmpdir=/work/kms173/tmp -jar ~/bin/picard.jar CleanSam \
 VALIDATION_STRINGENCY=LENIENT INPUT=$homedir/data/bam/${prefix}_sorted.bam OUTPUT=$homedir/data/bam/$prefix.clean.bam 2> $log/$prefix.cleansam.log
 
#Sort the file by the position
~/bin/jdk1.8.0_101/bin/java -Djava.io.tmpdir=/work/kms173/tmp -jar ~/bin/picard.jar SortSam \
VALIDATION_STRINGENCY=LENIENT INPUT=$homedir/data/bam/$prefix.clean.bam OUTPUT=$homedir/data/bam/$prefix.sort.bam SORT_ORDER=coordinate 2> $log/$prefix.sortsam.log

#Add read groups to data
~/bin/jdk1.8.0_101/bin/java -Djava.io.tmpdir=/work/kms173/tmp -jar ~/bin/picard.jar AddOrReplaceReadGroups \
I=$homedir/data/bam/$prefix.sort.bam O=$homedir/data/bam/$prefix.sortrg.bam SORT_ORDER=coordinate RGID=$prefix RGLB=$project RGPL=ILLUMINA RGPU=$project RGSM=$prefix CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT 2> $log/$prefix.addRG.log

rm $homedir/data/bam/$prefix.clean.bam
rm $homedir/data/bam/$prefix.sort.bam
rm $homedir/data/bam/${prefix}_sorted.bam
rm $homedir/data/bam/${prefix}_sorted.bai
