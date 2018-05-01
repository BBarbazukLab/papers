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

## (1) run on data where 1st 10 bases removed

DATAIN=/scratch/lfs/mcintyre/trago/outfiles/trimmomatic
OUTDIR=/scratch/lfs/mcintyre/trago/trago_output
OUT=table_trago_trimmomatic_cnts.txt

cd $DATAIN
for FILE in *.fq
do 
	CNT=`grep "@HWI" $FILE | wc -l`
	printf "$FILE,$CNT\n" >> $OUTDIR/$OUT.csv
done 

## (2) run on data with 10 bases kept

DATAIN2=/scratch/lfs/mcintyre/trago/outfiles/trimmomatic_10
OUTDIR=/scratch/lfs/mcintyre/trago/trago_output
OUT=table_trago_trimmomatic_10_cnts.txt

cd $DATAIN2
for FILE in *.fq
do 
  	CNT=`grep "@HWI" $FILE | wc -l`
        printf "$FILE,$CNT\n" >> $OUTDIR/$OUT.csv
done

