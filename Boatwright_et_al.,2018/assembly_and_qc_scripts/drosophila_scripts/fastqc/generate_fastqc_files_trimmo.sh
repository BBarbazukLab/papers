#!/bin/bash 

#generate fastqc files

PROJ=/bio/mcintyre/trago
DATAIN=$PROJ/outfiles/trimmomatic
OUT=$PROJ/outfiles/fastqc

cd $DATAIN

## output a list of files and step through that list.

ls > $PROJ/list_fq_files.txt 

for FILE in {1..48} #there are 48 files in the list
do

DESIGN=$(cat $PROJ/list_fq_files.txt | head -n $FILE | tail -n 1)  #steps through list_fq_files.txt by row

module load fastqc

fastqc $DESIGN --outdir=$OUT

done
