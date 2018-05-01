#!/bin/bash
#SBATCH --account=barbazuk
#SBATCH --qos=barbazuk
#SBATCH --mail-user=lboat@ufl.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/ufrc/barbazuk/lboat/Old_World_New_BED/scripts/PBS_LOGS/sam_compare.out
#SBATCH -e /ufrc/barbazuk/lboat/Old_World_New_BED/scripts/PBS_LOGS/sam_compare.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH -p bigmem
#SBATCH --mem-per-cpu=160G
#SBATCH --time=45:00:00

module load python/2.7.3

#### Set Directories
PROJ=/ufrc/barbazuk/lboat/Old_World_New_BED
SCRIPTS=/ufrc/barbazuk/lboat/Old_World_New_BED/scripts
TMPDIR=/ufrc/barbazuk/lboat/Old_World_New_BED/temp_dir

OUTPUT=$PROJ/trago_output/ase_counts/parents
if [ ! -e $OUTPUT ]; then mkdir $OUTPUT; fi

#### Create LOG directory
LOGS=$OUTPUT/logs
if [ ! -e $LOGS ]; then mkdir $LOGS; fi

## TDU reads to TPU and TPO for Tm
for i in 1 2 3
do
## Doing Tdu -- TDU and TPO alignments for comparing
SAMA=$PROJ/trago_output/filtered_parent_sam_files/Tdu_${i}_unq_2_TDU_4_Tm_commonID.sam
SAMB=$PROJ/trago_output/filtered_parent_sam_files/Tdu_${i}_unq_2_TPO_4_Tm_commonID.sam

### Concatenate reads for use in sam-compare
cat ${PROJ}/original_data/Tdu${i}_*.fastq | sed 's/ /_/g' > $TMPDIR/Tdu${i}.fastq

python $SCRIPTS/sam_compare.py \
-d \
-l 100 \
-f $PROJ/references/TDU_tpo_bed_for_sam_compare.bed \
-q $TMPDIR/Tdu${i}.fastq \
-A $SAMA \
-B $SAMB \
-c $OUTPUT/ase_counts_Tdu_${i}_2_Tdu_Tpo.csv \
-t $OUTPUT/ase_totals_Tdu_${i}_2_Tdu_Tpo.csv \
-g $LOGS/ase_counts_Tdu_${i}_2_Tdu_Tpo.log
done


## TDU reads to TPU and TPR for Tms
for i in 1 2 3
do

## Doing Tdu -- TDU and TPR alignments for comparing
SAMA=$PROJ/trago_output/filtered_parent_sam_files/Tdu_${i}_unq_2_TDU_4_Tms_commonID.sam
SAMB=$PROJ/trago_output/filtered_parent_sam_files/Tdu_${i}_unq_2_TPR_4_Tms_commonID.sam

### Concatenate reads for use in sam-compare
cat ${PROJ}/original_data/Tdu${i}_*.fastq | sed 's/ /_/g' > $TMPDIR/Tdu${i}.fastq

python $SCRIPTS/sam_compare.py \
-d \
-l 100 \
-f $PROJ/references/TDU_tpr_bed_for_sam_compare.bed \
-q $TMPDIR/Tdu${i}.fastq \
-A $SAMA \
-B $SAMB \
-c $OUTPUT/ase_counts_Tdu_${i}_2_Tdu_Tpr.csv \
-t $OUTPUT/ase_totals_Tdu_${i}_2_Tdu_Tpr.csv \
-g $LOGS/ase_counts_Tdu_${i}_2_Tdu_Tpr.log
done


## TPO reads to TDU and TPO for Tm
for i in 1 2 3
do

## Doing Tpo -- TDU and TPO alignments for comparing
SAMA=$PROJ/trago_output/filtered_parent_sam_files/Tpo_${i}_unq_2_TDU_4_Tm_commonID.sam
SAMB=$PROJ/trago_output/filtered_parent_sam_files/Tpo_${i}_unq_2_TPO_4_Tm_commonID.sam

### Concatenate reads for use in sam-compare
cat ${PROJ}/original_data/Tpo${i}_*.fastq | sed 's/ /_/g' > $TMPDIR/Tpo${i}.fastq

python $SCRIPTS/sam_compare.py \
-d \
-l 100 \
-f $PROJ/references/TPO_tdu_bed_for_sam_compare.bed \
-q $TMPDIR/Tpo${i}.fastq \
-A $SAMA \
-B $SAMB \
-c $OUTPUT/ase_counts_Tpo_${i}_2_Tdu_Tpo.csv \
-t $OUTPUT/ase_totals_Tpo_${i}_2_Tdu_Tpo.csv \
-g $LOGS/ase_counts_Tpo_${i}_2_Tdu_Tpo.log
done

## TPR reads to TDU and TPR for Tms
for i in 1 2 3
do

## Doing Tpr -- TDU and TPR alignments for comparing
SAMA=$PROJ/trago_output/filtered_parent_sam_files/Tpr_${i}_unq_2_TDU_4_Tms_commonID.sam
SAMB=$PROJ/trago_output/filtered_parent_sam_files/Tpr_${i}_unq_2_TPR_4_Tms_commonID.sam

### Concatenate reads for use in sam-compare
cat ${PROJ}/original_data/Tpr${i}_*.fastq | sed 's/ /_/g' > $TMPDIR/Tpr${i}.fastq

python $SCRIPTS/sam_compare.py \
-d \
-l 100 \
-f $PROJ/references/TPR_tdu_bed_for_sam_compare.bed \
-q $TMPDIR/Tpr${i}.fastq \
-A $SAMA \
-B $SAMB \
-c $OUTPUT/ase_counts_Tpr_${i}_2_Tdu_Tpr.csv \
-t $OUTPUT/ase_totals_Tpr_${i}_2_Tdu_Tpr.csv \
-g $LOGS/ase_counts_Tpr_${i}_2_Tdu_Tpr.log
done
