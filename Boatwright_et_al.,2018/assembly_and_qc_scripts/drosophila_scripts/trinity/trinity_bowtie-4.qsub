#!/bin/bash
#PBS -N Drosophila
#PBS -M anazarian@ufl.edu
#PBS -m n
#PBS -r n
#PBS -o /bio/mcintyre/trago/outfiles/PBS_LOGS/trinity_bowtie
#PBS -o Drosophila_trinity_bowtie_4.log
#PBS -j oe
#PBS -l nodes=1:ppn=4
#PBS -l walltime=15:00:00
#PBS -l pmem=5gb
#PBS -q bio
#PBS -t 1

#### This specifies to use the directory I submitted the script from
cd $PBS_O_WORKDIR
date
module load trinity

####module load bowtie

#### Set Directories
WORK=/bio/mcintyre
PROJ=trago
INDIR=$WORK/$PROJ/outfiles/trimmomatic_concatfiles
OUTDIR=$WORK/$PROJ/outfiles
if [ ! -e $OUTDIR ]; then mkdir -p $OUTDIR; fi

#### Because I am using an Array I am pulling LINE:MV:REP from an external CSV with all possible combinations

DESIGN_FILE=$WORK/$PROJ/dros_biorep_list_2.csv
DESIGN=$(sed -n "${PBS_ARRAYID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"
LINE=${ARRAY[0]}
MV=${ARRAY[1]}
REP=${ARRAY[2]}
LANE1=${ARRAY[3]}
LANE2=`echo $LANE1 | sed 's/\.1/\.2/g'`

#### Create LOG directory and start log
LOGS=$OUTDIR/LOGS/trinity_bowtie 		#script log information
if [ ! -e $LOGS ]
then
    mkdir -p $LOGS
fi

MYLOG=$LOGS/${LINE}_${MV}_${REP}_${LANE1}_${LANE2}_Drosophila_trinity_bowtie_4.log
printf "`date` $LINE $MV $REP $LANE1 $LANE2 SGE_ID:$PBS_ARRAYID HOSTNAME:$HOSTNAME \n\n" > $MYLOG

#####TRINLOG=$LOGS/${LINE}_${MV}_${REP}_${LANE1}_${LANE2}_Drosophila_trinity.log

#### Start Bowtie script
READS1=$INDIR/${LINE}_${MV}_${REP}_${LANE1}_trimmomatic_paired.fq
READS2=$INDIR/${LINE}_${MV}_${REP}_${LANE2}_trimmomatic_paired.fq
TARGET=$WORK/$PROJ/outfiles/trinity/Trinity.fasta
OUT1=$OUTDIR/trinity_bowtie_4
if [ ! -e $OUT1 ]
then
    mkdir -p $OUT1
fi

#### RUN Bowtie on trinity outfile and Paired-end Files

/apps/trinity/r20130225/util/alignReads.pl --left $READS1 --right $READS2 --seqType fq --target $TARGET --aligner bowtie --retain_SAM_file --retain_intermediate_files --output $OUT1 --prep_rsem -- -p 4 --strata --tryhard --best -n 3 -m 1 -l 28


2>>$MYLOG
