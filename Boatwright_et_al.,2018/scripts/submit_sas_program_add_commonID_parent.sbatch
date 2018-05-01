#!/bin/bash
#SBATCH --account=barbazuk
#SBATCH --qos=barbazuk
#SBATCH --mail-user=lboat@ufl.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/ufrc/barbazuk/lboat/Old_World_New_BED/scripts/PBS_LOGS/add_commonID.out
#SBATCH -e /ufrc/barbazuk/lboat/Old_World_New_BED/scripts/PBS_LOGS/add_commonID.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=45:00:00
module load sas/9.4

PROJ=/ufrc/barbazuk/lboat/Old_World_New_BED
LOGS=${PROJ}/sas_programs/logs
OUTPUT=${PROJ}/sas_programs/output
TMPDIR=/ufrc/barbazuk/lboat/Old_World_New_BED/sas_temp


## add commonID to parent sam files
sas -log ${LOGS}/add_commonID_to_parent_sam_files_HPC.log \
-print ${OUTPUT}/add_commonID_to_parent_sam_files_HPC.prt \
-work $TMPDIR \
-sysin ${PROJ}/sas_programs/add_commonID_to_parent_sam_files_HPC.sas
