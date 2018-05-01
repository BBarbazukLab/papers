#!/bin/bash

module load sas/9.4

PROJ=/ufrc/barbazuk/lboat/Old_World_New_BED
LOGS=$PROJ/sas_data/logs
OUTPUT=$PROJ/sas_data/output

## fix bed file coordinates
sas -log $LOGS/import_bayesian_results_hybrid.log \
    -print $OUTPUT/import_bayesian_results_hybrid.prt \
    -work $PROJ/sas_downstream_scripts/tmpdir \
    -sysin $PROJ/sas_downstream_scripts/import_bayesian_results_hybrid.sas

