#!/bin/bash
#PBS -N Drosophila
#PBS -M anazarian@ufl.edu
#PBS -m n
#PBS -r n
#PBS -o /bio/mcintyre/trago/outfiles/PBS_LOGS/trimmomatic_concatfiles
#PBS -o Drosophila_trimmomatic_concatfiles.log
#PBS -j oe
#PBS -l nodes=1:ppn=1
#PBS -l walltime=25:00:00
#PBS -l pmem=5gb
#PBS -q bio
#PBS -t 1%1

#### This specifies to use the directory I submitted the script from
cd $PBS_O_WORKDIR
date
module load trimmomatic

#### Set Directories
WORK=/bio/mcintyre
PROJ=trago
INDIR=$WORK/$PROJ/outfiles/cutadapt_concatfiles
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
LOGS=$OUTDIR/LOGS/trimmomatic_concatfiles 		#script log information
if [ ! -e $LOGS ]
then
    mkdir -p $LOGS
fi

####MYLOG=$LOGS/${LINE}_${MV}_${REP}_${LANE}_trimmomatic.log
####printf "`date` $LINE $MV $REP SGE_ID:$PBS_ARRAYID HOSTNAME:$HOSTNAME \n\n" > $MYLOG

TRIMOLOG=$LOGS/${LINE}_${MV}_${REP}_${LANE1}_${LANE2}_Drosophila_trimmomatic.log

#### Start Trimmomatic
READS1=$INDIR/${LINE}_${MV}_${REP}_${LANE1}_adapter_trimmed.fq
READS2=$INDIR/${LINE}_${MV}_${REP}_${LANE2}_adapter_trimmed.fq
OUT1=$OUTDIR/trimmomatic_concatfiles
if [ ! -e $OUT1 ]
then
    mkdir -p $OUT1
fi

TRIM1=$OUT1/${LINE}_${MV}_${REP}_${LANE1}_trimmomatic_paired.fq
TRIM2=$OUT1/${LINE}_${MV}_${REP}_${LANE1}_trimmomatic_unpaired.fq
TRIM3=$OUT1/${LINE}_${MV}_${REP}_${LANE2}_trimmomatic_paired.fq
TRIM4=$OUT1/${LINE}_${MV}_${REP}_${LANE2}_trimmomatic_unpaired.fq

### Start Lane Section in Log
####printf "$LANE\n" >> $MYLOG

#### RUN trimmomatic on Paired end Files
java -classpath $HPC_TRIMMOMATIC_LIB/trimmomatic.jar org.usadellab.trimmomatic.TrimmomaticPE -phred64 -threads 3 -trimlog $TRIMOLOG $READS1 $READS2 $TRIM1 $TRIM2 $TRIM3 $TRIM4 HEADCROP:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:60	

####2>>$MYLOG
