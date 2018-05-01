
module load sas/9.4

PROJ=/ufrc/barbazuk/lboat/Old_World_New_BED
LOGS=$PROJ/sas_data/logs
OUTPUT=$PROJ/sas_data/output
TMPDIR=$PROJ/temp_dir/

## fix bed file coordinates
sas -log $LOGS/import_ase_counts_hybrid.log \
    -print $OUTPUT/import_ase_counts_hybrid.prt \
    -work $TMPDIR \
    -sysin $PROJ/sas_downstream_scripts/import_ase_counts_hybrid.sas

