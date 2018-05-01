#!/bin/bash
#SBATCH --account=barbazuk
#SBATCH --qos=barbazuk-b
#SBATCH --mail-user=lboat@ufl.edu
#SBATCH --job-name=sam_hybrid
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH -o /ufrc/barbazuk/lboat/Old_World_New_BED/scripts/PBS_LOGS/hybrid_compare_%j.out
#SBATCH -e /ufrc/barbazuk/lboat/Old_World_New_BED/scripts/PBS_LOGS/hybrid_compare_%j.err
#SBATCH -t 0-30:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50G
#SBATCH --partition=bigmem

module load python/2.7.3 

#### Set Directories
PROJ=/ufrc/barbazuk/lboat/Old_World_New_BED
SCRIPTS=$PROJ/scripts
ORIG=$PROJ/original_data/trago-tetraploids
TMPDIR=$PROJ/temp_dir

OUTPUT=$PROJ/trago_output/ase_counts/hybrids
if [ ! -e $OUTPUT ]; then mkdir $OUTPUT; fi

#### Create LOG directory
LOGS=$OUTPUT/logs
if [ ! -e $LOGS ]; then mkdir $LOGS; fi

## Cast reads to Lam and Croc 
for i in {1..3} 
do
    ## Doing CAST -- Lam and Croc alignments for comparing
    SAMA=$PROJ/trago_output/filtered_hybrid_sam_files/Tm_${i}_unq_2_TDU_commonID.sam
    SAMB=$PROJ/trago_output/filtered_hybrid_sam_files/Tm_${i}_unq_2_TPO_commonID.sam

    ### Concatenate reads for use in sam-compare
    cat $ORIG/Sample_Tm${i}/*.fastq | sed 's/ /_/g' > $TMPDIR/Tm_${i}/Tm$i.fastq

    python $SCRIPTS/sam_compare.py \
        -d \
	-l 100 \
        -f $PROJ/references/TDU_tpo_bed_for_sam_compare.bed \
        -q $TMPDIR/Tm_${i}/Tm${i}.fastq \
        -A $SAMA \
        -B $SAMB \
        -c $OUTPUT/ase_counts_Tm_${i}_2_tdu_tpo.csv \
        -t $OUTPUT/ase_totals_Tm_${i}_2_tdu_tpo.csv \
        -g $LOGS/ase_counts_Tm_${i}_2_tdu_tpo.log
done


for i in {1..3}
do
    ## Doing CAST -- Lam and Croc alignments for comparing
    SAMA=$PROJ/trago_output/filtered_hybrid_sam_files/Tms_${i}_unq_2_TDU_commonID.sam
    SAMB=$PROJ/trago_output/filtered_hybrid_sam_files/Tms_${i}_unq_2_TPR_commonID.sam

    ### Concatenate reads for use in sam-compare
    cat $ORIG/Sample_Tms${i}/*.fastq | sed 's/ /_/g' > $TMPDIR/Tms_${i}/Tms$i.fastq

    python $SCRIPTS/sam_compare.py \
        -d \
	-l 100 \
        -f $PROJ/references/TDU_tpr_bed_for_sam_compare.bed \
        -q $TMPDIR/Tms_${i}/Tms${i}.fastq \
        -A $SAMA \
        -B $SAMB \
        -c $OUTPUT/ase_counts_Tms_${i}_2_tdu_tpr.csv \
        -t $OUTPUT/ase_totals_Tms_${i}_2_tdu_tpr.csv \
        -g $LOGS/ase_counts_Tms_${i}_2_tdu_tpr.log
done

