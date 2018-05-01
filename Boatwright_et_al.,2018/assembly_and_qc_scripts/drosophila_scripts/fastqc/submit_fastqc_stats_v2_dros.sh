#!/bin/bash
#PBS -M anazarian@ufl.edu
#PBS -m n
#PBS -r n
#PBS -q bio
#PBS -l walltime=1:00:00
#PBS -l nodes=1:ppn=1
#PBS -l pmem=1gb
#PBS -j oe
#PBS -o /project/mcintyre/trago


PROJ=/project/mcintyre/trago/dros_out
DATAIN=$PROJ/fastqc
#####OUTDIR=$PROJ/fastqc_summary/files
OUTDIR=$PROJ/summary-fq

cd $DATAIN

for DIR in $( find ./ -maxdepth 1 -type d | cut -f2 -d'/')
do 

	cd $DATAIN/$DIR

	NAME=$(basename "$DIR" _fastqc)

	    perl /project/mcintyre/dandelion_2012_HiSeq/scripts/fastqc_stats.pl $NAME $OUTDIR

done

