#!/bin/bash
#SBATCH --account=barbazuk
#SBATCH --qos=barbazuk-b
#SBATCH --mail-user=lboat@ufl.edu
#SBATCH --job-name=filter_sam
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH -o /ufrc/barbazuk/lboat/Old_World_New_BED/scripts/PBS_LOGS/filter_sam_cast_2_parents_%a.out
#SBATCH -e /ufrc/barbazuk/lboat/Old_World_New_BED/scripts/PBS_LOGS/filter_sam_cast_2_parents_%a.err
#SBATCH -t 0-24:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --mem-per-cpu=5G 

module load perl/5.16.0

#### Set Directories
PROJ=/ufrc/barbazuk/lboat/Old_World_New_BED
SCRIPTS=$PROJ/scripts
REF=$PROJ/references
SAM=$PROJ/trago_output/aln_hybrids_to_consensed

OUTPUT=$PROJ/trago_output/filtered_hybrid_sam_files 
if [ ! -e $OUTPUT ]; then mkdir $OUTPUT; fi

#### Create LOG directory
LOGS=$OUTPUT/logs
if [ ! -e $LOGS ]; then mkdir $LOGS; fi


# 3 reps for Tm
for REP in {1..3}
do

# filter sam file of Lam reads to Lam with Lam bed file
    perl $SCRIPTS/sam-filter-by-bed.pl \
        -b $REF/TDU-tpo_overlaps_WRT_orthologs.bed \
        -s $SAM/Tm_${REP}_uniq_2_TDU.sam \
        >$OUTPUT/Tm_${REP}_uniq_2_TDU_filter_for_Tm.sam \
        2> $LOGS/Tm_${REP}_uniq_2_TDU_filter_for_Tm.log

# filter sam file of Lam reads to Croc with Croc bed file
    perl $SCRIPTS/sam-filter-by-bed.pl \
        -b $REF/TPO-tdu_overlaps_WRT_orthologs.bed \
        -s $SAM/Tm_${REP}_uniq_2_TPO.sam \
        >$OUTPUT/Tm_${REP}_uniq_2_TPO_filter_for_Tm.sam \
        2> $LOGS/Tm_${REP}_uniq_2_TPO_filter_for_Tm.log

done

# 3 reps for Tms
for REP in {1..3}
do 

# filter sam file of Lam reads to Lam with Lam bed file
    perl $SCRIPTS/sam-filter-by-bed.pl \
	-b $REF/TDU-tpr_overlaps_WRT_orthologs.bed \
	-s $SAM/Tms_${REP}_uniq_2_TDU.sam \
	>$OUTPUT/Tms_${REP}_uniq_2_TDU_filter_for_Tms.sam \
	2> $LOGS/Tms_${REP}_uniq_2_TDU_filter_for_Tms.log

# filter sam file of Lam reads to Croc with Croc bed file
    perl $SCRIPTS/sam-filter-by-bed.pl \
        -b $REF/TPR-tdu_overlaps_WRT_orthologs.bed \
        -s $SAM/Tms_${REP}_uniq_2_TPR.sam \
        >$OUTPUT/Tms_${REP}_uniq_2_TPR_filter_for_Tms.sam \
        2> $LOGS/Tms_${REP}_uniq_2_TPR_filter_for_Tms.log

done



