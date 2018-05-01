#!/bin/bash
#SBATCH --mail-user=lboat@ufl.edu
#SBATCH --mail-type=ALL
#SBATCH -o /ufrc/barbazuk/lboat/Old_World_New_BED/scripts/PBS_LOGS/aln_Tms_2_lucas_TPR_consensed_contigs_%a.out
#SBATCH -e /ufrc/barbazuk/lboat/Old_World_New_BED/scripts/PBS_LOGS/aln_Tms_2_lucas_TPR_consensed_contigs_%a.err
#SBATCH -t 0-48:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 4

#SBATCH --mem-per-cpu=60gb
#SBATCH -p bigmem
#SBATCH --array=4-6
     ## 1-3 are Tm and 4-6 are Tms
     ## Tm = TDU x TPO
     ## Tms = TDU x TPR

#### align Tms reads to TPR consensed 

module load bowtie/0.12.9
module load last/531
module load python/2.7.3

#### Store number of processors used -- i.e. number of files to split into to run LAST
    NUMPROCS=4

## Set directories
    PROJ=/ufrc/barbazuk/lboat/Old_World_New_BED
    ORIG=$PROJ/original_data

    OUTPUT=$PROJ/trago_output/aln_hybrids_to_consensed
    if [ ! -e $OUTPUT ]; then mkdir -p $OUTPUT; fi
    ALN_LOGS=$OUTPUT/aln_logs
    if [ ! -e $ALN_LOGS ]; then mkdir -p $ALN_LOGS; fi

## Pull info from design file
    DESIGN_FILE=$PROJ/design_files/trago_hybrids_with_reps.txt
    DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

    IFS=',' read -ra ARRAY <<< "$DESIGN"

    SP=${ARRAY[0]}
    REP=${ARRAY[1]}

    TMPDIR=/ufrc/barbazuk/lboat/Old_World_New_BED/temp_dir/TPR_ref/${SP}_${REP}

## Reference
    REF=$PROJ/references/TPR_consensed_contigs_500-15k
    LASTREF=$PROJ/references/TPR_consensed_contigs_500-15k_last531

#### Create LOG directory and start log
    LOGS=$OUTPUT/job_logs              #script log information
    ALN_ERROR_LOG=$LOGS/size_errors.txt
    if [ ! -d $LOGS ]; then mkdir -p $LOGS; fi

    MYLOG=$LOGS/${SP}_${REP}_TPR_consensus.log
    printf "`date` $SP SGE_ID:$SLURM_ARRAY_TASK_ID HOSTNAME:$HOSTNAME \n\n" > $MYLOG

    ALNLOGS=$OUTPUT/aln_logs       
    if [ ! -d $ALNLOGS ]; then mkdir -p $ALNLOGS; fi


#### for each species, cat together R1 and R2 (treating as single end) and lanes (reps separate)
    cat ${PROJ}/original_data/trago-tetraploids/Sample_${SP}${REP}/${SP}${REP}_*.fastq | sed 's/ /_/g' > $TMPDIR/${SP}_${REP}.fastq

    READS=$TMPDIR/${SP}_${REP}.fastq

#### Determine my quality FASTQ score type
    qual=`python /ufrc/barbazuk/lboat/Old_World_New_BED/scripts/identify_quality.py -i $READS`
    if [ $qual == 'phred64' ];
    then
        # set to old illumina quality scores phred64/solexa 1.3
        btqual='--phred64-quals'
        lastqual='3'
    else
        # change to sanger format which is what all new illumina data is
        btqual='--phred33-quals'
        lastqual='1'
        fi

    cd $TMPDIR

#### run bowtie uniq
    echo "START Single End bowtie with 3 mm [`date`]" >>$MYLOG
    bowtie -S \
           $btqual \
           --best \
           --tryhard \
           --strata \
           --chunkmbs 1024 \
           --sam-no-qname-trunc \
           -p $NUMPROCS \
           -m 1 \
           -v 3 \
           --un ${SP}_${REP}_unaln_bt.fq \
           --max ${SP}_${REP}_ambig_bt.fq \
           $REF \
           $READS \
           2>$ALNLOGS/${SP}_${REP}_2_TPR_bowtie.log \
           | perl -ne 'if(/^@/){next} @s = split("\t"); if($s[1] != 4) {print}' \
           > ${SP}_${REP}_aln_bt.sam \
           2>>$MYLOG
        echo "FINISH Single End bowtie [`date`]" >>$MYLOG

        if [ ! -e ${SP}_${REP}_unaln_bt.fq ]
        then
            echo "WARNING: There were no unaligned reads" >> $MYLOG
        fi

        cat ${SP}_${REP}_unaln_bt.fq  ${SP}_${REP}_ambig_bt.fq > ${SP}_${REP}.fq
        READS=${SP}_${REP}.fq

#### LAST all

    # The following Python script turns READS into NUMPROCS separate files, going from 0-split.fq to (NUMPROCS-1)-split.fq. 
    # READS is the location where the output files are placed

        python /ufrc/barbazuk/lboat/Old_World_New_BED/scripts/SeparateLast.py $READS $NUMPROCS

    # START LAST for each of the NUMPROCS new -split files
        echo "START Single End LAST [`date`]" >>$MYLOG
        for ((I=0; I<$NUMPROCS; I++))
        do
                READS=${I}-split.fq

                lastal \
                    -l 25 \
                    -Q $lastqual \
                    $LASTREF \
                    $READS \
                    2>>$MYLOG \
                    >${I}_last.maf &
        done
        wait

        echo "FINISH Single End LAST [`date`]" >>$MYLOG

        rm *-split.fq

    # Convert MAF to SAM
        echo "START Converting from MAF to SAM [`date`]" >>$MYLOG
        for ((I=0; I<$NUMPROCS; I++))
        do
            python $PROJ/scripts/maf-convert.py sam ${I}_last.maf > ${I}_last.sam 2>>$MYLOG &
        done
        wait

        rm *_last.maf

        echo "FINISH Converting to SAM [`date`]" >>$MYLOG

        for ((I=0; I<$NUMPROCS; I++))
        do
            cat ${I}_last.sam >> ${SP}_${REP}_last.sam
            rm ${I}_last.sam
        done

    # PARSE LAST FILES
        echo "START Parsing LAST SAM file into unique and ambiguous [`date`]" >>$MYLOG
        perl /ufrc/barbazuk/lboat/Old_World_New_BED/scripts/parse_last_sam_v2.pl \
            ${SP}_${REP}_last.sam \
            ${SP}_${REP}_ambig_last.sam \
            ${SP}_${REP}_uniq_last.sam \
            > $ALNLOGS/${SP}_${REP}_2_TPR_LAST.log
            rm ${SP}_${REP}_last.sam 

        echo "FINISH Parsing LAST SAM file into unique and ambiguous [`date`]" >>$MYLOG

#### File Check
    cd $TMPDIR
    echo "FILE SIZE CHECK">>$MYLOG
    du -sh *.sam >>$MYLOG

#### Combine all Sam files
    echo "START Combine SAM files">>$MYLOG

    cat ${SP}_${REP}_aln_bt.sam ${SP}_${REP}_uniq_last.sam > $OUTPUT/${SP}_${REP}_uniq_2_TPR.sam 2>>$MYLOG
    mv ${SP}_${REP}_ambig_last.sam $OUTPUT/${SP}_${REP}_ambig_last_2_TPR.sam 2>>$MYLOG
    
echo "Script complete, [`date`]" >>$MYLOG


