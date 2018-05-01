#!/bin/bash

module load sas/9.4

PROJ=/ufrc/barbazuk/lboat/Old_World_New_BED
LOGS=$PROJ/sas_data/logs
OUTPUT=$PROJ/sas_data/output

## fix bed file coordinates
sas -log $LOGS/bayesian_flag_sig_hybrid.log \
    -print $OUTPUT/bayesian_flag_sig_hybird.prt \
    -work $PROJ/sas_downstream_scripts/tmpdir \
    -sysin $PROJ/sas_downstream_scripts/bayesian_flag_sig_results_hybrid_pdf.sas

