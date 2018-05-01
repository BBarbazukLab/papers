#!/bin/bash
#PBS -N Drosophila
#PBS -M anazarian@ufl.edu
#PBS -m n
#PBS -r n
#PBS -o /bio/mcintyre/trago/dros_outfiles/PBS_LOGS/bowtie
#PBS -o Drosophila_bowtie.log
#PBS -j oe
#PBS -l nodes=1:ppn=4
#PBS -l walltime=10:00:00
#PBS -l pmem=5gb
#PBS -q bio
#PBS -t 1


#### This specifies to use the directory I submitted the script from
cd $PBS_O_WORKDIR
date
module load bowtie

####module load bowtie

#### Set Directories
WORK=/bio/mcintyre
PROJ=trago
INDIR=$WORK/$PROJ/dros_outfiles/trimmomatic_concatfiles
OUTDIR=$WORK/$PROJ/dros_outfiles
if [ ! -e $OUTDIR ]; then mkdir -p $OUTDIR; fi

#### Because I am using an Array I am pulling LINE:MV:REP from an external CSV with all possible combinations

DESIGN_FILE=$WORK/$PROJ/dros_outfiles/dros_biorep_list_2.csv
DESIGN=$(sed -n "${PBS_ARRAYID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"
LINE=${ARRAY[0]}
MV=${ARRAY[1]}
REP=${ARRAY[2]}
LANE1=${ARRAY[3]}
LANE2=`echo $LANE1 | sed 's/\.1/\.2/g'`

#### Create LOG directory and start log
LOGS=$OUTDIR/LOGS/bowtie 		#script log information
if [ ! -e $LOGS ]
then
    mkdir -p $LOGS
fi

MYLOG=$LOGS/${LINE}_${MV}_${REP}_${LANE1}_${LANE2}_Drosophila_bowtie.log
printf "`date` $LINE $MV $REP $LANE1 $LANE2 SGE_ID:$PBS_ARRAYID HOSTNAME:$HOSTNAME \n\n" > $MYLOG

#### Start Bowtie script
READS1=$INDIR/${LINE}_${MV}_${REP}_${LANE1}_trimmomatic_paired.fq
READS2=$INDIR/${LINE}_${MV}_${REP}_${LANE2}_trimmomatic_paired.fq
TARGET=$WORK/$PROJ/dros_outfiles/references/dmel-all-transcript-r5.30.fasta
OUT1=$OUTDIR/bowtie
if [ ! -e $OUT1 ]
then
    mkdir -p $OUT1
fi

REF=$OUT1/bowtie_REF
if [ ! -e $REF ]
then
    mkdir -p $REF
fi


#### RUN Bowtie on trinity outfile and Paired-end Files

bowtie-build $TARGET $REF/dros_ref
bowtie --chunkmbs 1040 --sam -n 3 -e 70 -l 29 -m 1 -p 4 --best --strata --tryhard --un $OUT1/unaln_dros.fq --max $OUT1/ambig_dros.fq $REF/dros_ref -1 $READS1 -2 $READS2 | perl -ne 'if(/^@/){next} @s = split; if($s[1] != 77 && $s[1] != 141) {print}' > $OUT1/ALN_dros.sam

2>>$MYLOG
