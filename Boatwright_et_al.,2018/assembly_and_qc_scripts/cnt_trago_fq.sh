#!/bin/bash
#PBS -M ammorse@ufl.edu
#PBS -m n
#PBS -r n
#PBS -l walltime=1:00:00
#PBS -l nodes=1:ppn=1
#PBS -l pmem=1gb
#PBS -o /scratch/lfs/mcintyre/trago/scripts/PBS_LOGS/
#PBS -j oe
#PBS -q bio

## Set directories
ORIG=/scratch/lfs/mcintyre/trago/trago_data
OUTDIR=/scratch/lfs/mcintyre/trago/trago_output
OUT=table_trago_fq_cnts.txt

cd $ORIG
for FILE in *.fastq
do 
	CNT=`grep "@HWI" $FILE | wc -l`
	printf "$FILE,$CNT\n" >> $OUTDIR/$OUT.csv
done 


