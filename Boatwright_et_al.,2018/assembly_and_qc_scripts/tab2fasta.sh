#!/bin/bash
#PBS -M ammorse@ufl.edu
#PBS -m n
#PBS -q bio
#PBS -r n
#PBS -j oe
#PBS -o /scratch/lfs/mcintyre/cegs/scripts/PBS_LOGS/
#PBS -l walltime=1:00:00
#PBS -l nodes=1:ppn=1
#PBS -l pmem=3Gb


PROJ=/scratch/lfs/mcintyre/trago
DATAIN=$PROJ/trago_output

## convert tab to fasta format
	## 1st column is identifier
	## 2nd column in sequence

cd $DATAIN
for FILE in *table_uniq.txt 
do 
	cat -n $FILE | awk '{print ">"$1"_"$2"\n"$3}'  >$FILE.fa
done 

#cat -n Tdu_R1_table_uniq.txt | awk '{print ">"$1"_"$2"\n"$3}'  >Tdu_R1_uniq.fa
#cat -n Tdu_R2_table_uniq.txt | awk '{print ">"$1"_"$2"\n"$3}'  >Tdu_R2_uniq.fa

#cat -n Tpo_R1_table_uniq.txt | awk '{print ">"$1"_"$2"\n"$3}'  >Tpo_R1_uniq.fa
#cat -n Tpr_R2_table_uniq.txt | awk '{print ">"$1"_"$2"\n"$3}'  >Tpr_R2_uniq.fa

#cat -n Tpr_R1_table_uniq.txt | awk '{print ">"$1"_"$2"\n"$3}'  >Tpr_R1_uniq.fa
#cat -n Tpr_R2_table_uniq.txt | awk '{print ">"$1"_"$2"\n"$3}'  >Tpr_R2_uniq.fa


