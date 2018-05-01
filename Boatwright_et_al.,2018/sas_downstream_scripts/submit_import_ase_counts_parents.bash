module load sas

PROJ=/ufrc/barbazuk/lboat/Old_World_New_BED
LOGS=$PROJ/sas_data/logs
OUTPUT=$PROJ/sas_data/output
TMPDIR=$PROJ/sas_data/temp/

## fix bed file coordinates
sas -log $LOGS/import_ase_counts_parents.log \
    -print $OUTPUT/import_ase_counts_parents.prt \
    -work $TMPDIR \
    -sysin $PROJ/sas_downstream_scripts/import_ase_counts_parents.sas

