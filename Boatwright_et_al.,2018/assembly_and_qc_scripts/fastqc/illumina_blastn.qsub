#!/bin/bash
#PBS -N Ilumina-Blast
#PBS -M polvadore@ufl.edu
#PBS -m n
#PBS -r n
#PBS -o /bio/mcintyre/trago/scripts/fastqc/PBS_LOGS/blast
#PBS -o Illumina.log
#PBS -j oe
#PBS -l nodes=1:ppn=2
#PBS -l walltime=03:00:00
#PBS -l pmem=5gb
#PBS -q bio
#PBS -t 1

#### This specifies to use the directory I submitted the script from
cd $PBS_O_WORKDIR
date
module load ncbi_blast

#### Set Directories
WORK=/bio/mcintyre
PROJ=trago
INDIR=$WORK/$PROJ/outfiles/illumina_blast/fasta
OUTDIR=$WORK/$PROJ/outfiles/illumina_blast
if [ ! -e $OUTDIR ]; then mkdir -p $OUTDIR; fi


#### Because I am using an Array I am pulling REP:BARCODE:LANE:READTYPE1:READTYPE2 from an external CSV with all possible combinations

DESIGN_FILE=$WORK/$PROJ/trago_biorep_list.csv
DESIGN=$(sed -n "${PBS_ARRAYID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"
REP=${ARRAY[0]}
BARCODE=${ARRAY[1]}
LANE=${ARRAY[2]}
READTYPE=${ARRAY[3]}

#### Create LOG directory and start log
LOGS=$OUTDIR/LOGS/illumina_blast           #script log information
if [ ! -e $LOGS ]
then
    mkdir -p $LOGS
fi

MYLOG=$LOGS/${REP}_${BARCODE}_${LANE}_${READTYPE}_001.log
printf "`date` $REP $BARCODE $LANE $READTYPE SGE_ID:$PBS_ARRAYID HOSTNAME:$HOSTNAME \n\n" > $MYLOG

REF=$WORK/references/illumina_adapters/illumina_adapters_BLAST
READS=$INDIR/${REP}_${BARCODE}_${LANE}_${READTYPE}_001.fasta
OUTFILE=$OUTDIR/${REP}_${BARCODE}_${LANE}_${READTYPE}_001_BLAST.tsv

#### Start BLAST

blastn -db $REF -query $READS -outfmt 7 > $OUTFILE

2>>$MYLOG
