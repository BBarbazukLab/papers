module load sas

PROJ=/ufrc/barbazuk/lboat/T.lamottei_T.crocifolius/Complete_data_set/align_reads_for_samcompare
LOGS=$PROJ/sas_data/logs
OUTPUT=$PROJ/sas_data/output
TMPDIR=/ufrc/barbazuk/lboat/T.lamottei_T.crocifolius/Complete_data_set/align_reads_for_samcompare/temp_files/

## fix bed file coordinates
sas -log $LOGS/import_ase_counts_parents_L3.log \
    -print $OUTPUT/import_ase_counts_parents_L3.prt \
    -work $TMPDIR \
    -sysin $PROJ/sas_downstream_scripts/import_ase_counts_parents_Lam3.sas

