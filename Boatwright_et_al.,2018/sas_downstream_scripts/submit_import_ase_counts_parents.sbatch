#!/bin/bash 
#SBATCH --mail-user=lboat@ufl.edu
#SBATCH --job-name=sas
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o /ufrc/barbazuk/lboat/T.lamottei_T.crocifolius/Complete_data_set/align_reads_for_samcompare/PBS_LOGS/
#SBATCH -t 0-15:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --qos=barbazuk-b
#SBATCH --mem-per-cpu=1G

module load sas/9.4

PROJ=/ufrc/barbazuk/lboat/T.lamottei_T.crocifolius/Complete_data_set/align_reads_for_samcompare
LOGS=$PROJ/sas_data/logs
OUTPUT=$PROJ/sas_data/output
TMPDIR=/ufrc/barbazuk/lboat/T.lamottei_T.crocifolius/Complete_data_set/align_reads_for_samcompare/temp_files/

cd $SLURM_SUBMIT_DIR

## fix bed file coordinates
sas -log $LOGS/import_ase_counts_parents.log \
    -print $OUTPUT/import_ase_counts_parents.prt \
    -work $TMPDIR \
    -sysin $PROJ/sas_downstream_scripts/import_ase_counts_parents.sas

