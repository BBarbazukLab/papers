#!/bin/bash 

module load sas/9.4

PROJ=/ufrc/barbazuk/lboat/Old_World_New_BED
LOGS=$PROJ/sas_data/logs
OUTPUT=$PROJ/sas_data/output
TMPDIR=/ufrc/barbazuk/lboat/Old_World_New_BED/temp_dir/


## add commonID to parent sam files
sas -log $LOGS/add_commonID_to_hybrid_sam_files_HPC.log \
    -print $OUTPUT/add_commonID_to_hybrid_sam_files_HPC.prt \
    -work $TMPDIR \
    -sysin $PROJ/scripts/add_commonID_to_hybrid_sam_files_HPC.sas

