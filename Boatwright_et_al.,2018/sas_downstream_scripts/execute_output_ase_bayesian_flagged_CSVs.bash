#!/bin/bash

module load sas/9.4

PROJ=/ufrc/barbazuk/lboat/Old_World_New_BED
LOGS=$PROJ/sas_data/logs
OUTPUT=$PROJ/sas_data/output

## fix bed file coordinates
sas -log $LOGS/output_ase_bayesian_flagged_CSVs.log \
    -print $OUTPUT/output_ase_bayesian_flagged_CSVs.prt \
    -work $PROJ/sas_downstream_scripts/tmpdir \
    -sysin $PROJ/sas_downstream_scripts/output_ase_bayesian_flagged_CSVs.sas

