#!/bin/bash
#SBATCH --account=barbazuk
#SBATCH --qos=barbazuk
#SBATCH --mail-user=lboat@ufl.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/ufrc/barbazuk/lboat/Old_World_New_BED/scripts/PBS_LOGS/filter_NW_SAMs.out
#SBATCH -e /ufrc/barbazuk/lboat/Old_World_New_BED/scripts/PBS_LOGS/filter_NW_SAMs.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=5G
#SBATCH --time=48:00:00

module load perl

#### Set Directories
PROJ=/ufrc/barbazuk/lboat/Old_World_New_BED
SCRIPTS=/ufrc/barbazuk/lboat/Old_World_New_BED/scripts
REF=/ufrc/barbazuk/lboat/Old_World_New_BED/references
SAM=/ufrc/barbazuk/lboat/Old_World_New_BED/trago_output/aln_parents_to_consensed

OUTPUT=$PROJ/trago_output/filtered_parent_sam_files
if [ ! -e $OUTPUT ]; then mkdir $OUTPUT; fi

#### Create LOG directory
LOGS=$OUTPUT/logs
if [ ! -e $LOGS ]; then mkdir $LOGS; fi


for REP in 1 2 3
do

#### Filter sam file of TPO reads aligned to TDU with TDU_tpo bed file
perl $SCRIPTS/sam-filter-by-bed.pl \
-b $REF/TDU-tpo_overlaps_WRT_orthologs.bed \
-s $SAM/Tpo_${REP}_uniq_2_TDU.sam \
>$OUTPUT/Tpo_${REP}_uniq_2_TDU_filter_for_Tm.sam \
2> $LOGS/Tpo_${REP}_uniq_2_TDU_filter_for_Tm.log

#### filter sam file of TPO reads aligned to TPO with TPO_tdu bed file
perl $SCRIPTS/sam-filter-by-bed.pl \
-b $REF/TPO-tdu_overlaps_WRT_orthologs.bed \
-s $SAM/Tpo_${REP}_uniq_2_TPO.sam \
>$OUTPUT/Tpo_${REP}_uniq_2_TPO_filter_for_Tm.sam \
2> $LOGS/Tpo_${REP}_uniq_2_TPO_filter_for_Tm.log
done

for REP in 1 2 3
do

#### Filter sam file of TPR reads aligned to TDU with TDU_tpr bed file
perl $SCRIPTS/sam-filter-by-bed.pl \
-b $REF/TDU-tpr_overlaps_WRT_orthologs.bed \
-s $SAM/Tpr_${REP}_uniq_2_TDU.sam \
>$OUTPUT/Tpr_${REP}_uniq_2_TDU_filter_for_Tms.sam \
2> $LOGS/Tpr_${REP}_uniq_2_TDU_filter_for_Tms.log

#### filter sam file of TPR reads aligned to TPR with TPR_tdu bed file
perl $SCRIPTS/sam-filter-by-bed.pl \
-b $REF/TPR-tdu_overlaps_WRT_orthologs.bed \
-s $SAM/Tpr_${REP}_uniq_2_TPR.sam \
>$OUTPUT/Tpr_${REP}_uniq_2_TPR_filter_for_Tms.sam \
2> $LOGS/Tpr_${REP}_uniq_2_TPR_filter_for_Tms.log

done

for REP in 1 2 3
do

#### Filter sam file of TDU reads aligned to TDU with TDU_tpo bed file
perl $SCRIPTS/sam-filter-by-bed.pl \
-b $REF/TDU-tpo_overlaps_WRT_orthologs.bed \
-s $SAM/Tdu_${REP}_uniq_2_TDU.sam \
>$OUTPUT/Tdu_${REP}_uniq_2_TDU_filter_for_Tm.sam \
2> $LOGS/Tdu_${REP}_uniq_2_TDU_filter_for_Tm.log

#### filter sam file of TDU reads aligned to TPO with TPO_tdu bed file
perl $SCRIPTS/sam-filter-by-bed.pl \
-b $REF/TPO-tdu_overlaps_WRT_orthologs.bed \
-s $SAM/Tdu_${REP}_uniq_2_TPO.sam \
>$OUTPUT/Tdu_${REP}_uniq_2_TPO_filter_for_Tm.sam \
2> $LOGS/Tdu_${REP}_uniq_2_TPO_filter_for_Tm.log

#### filter sam file of TDU reads aligned to TDU with TDU_tpr bed file
perl $SCRIPTS/sam-filter-by-bed.pl \
-b $REF/TDU-tpr_overlaps_WRT_orthologs.bed \
-s $SAM/Tdu_${REP}_uniq_2_TDU.sam \
>$OUTPUT/Tdu_${REP}_uniq_2_TDU_filter_for_Tms.sam \
2> $LOGS/Tdu_${REP}_uniq_2_TDU_filter_for_Tms.log

#### filter sam file of TDU reads aligned to TPR with TPR_tdu bed file
perl $SCRIPTS/sam-filter-by-bed.pl \
-b $REF/TPR-tdu_overlaps_WRT_orthologs.bed \
-s $SAM/Tdu_${REP}_uniq_2_TPR.sam \
>$OUTPUT/Tdu_${REP}_uniq_2_TPR_filter_for_Tms.sam \
2> $LOGS/Tdu_${REP}_uniq_2_TPR_filter_for_Tms.log

done

