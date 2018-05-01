#!/bin/bash - 
#===============================================================================
# 
#   DESCRIPTION: Split the ASE table into chunks of rows for faster processing
#   in the Bayesian machine.
# 
#===============================================================================

# Load python
module load python/2.7.6

# Set Variables
PROJ=/ufrc/barbazuk/lboat/Old_World_New_BED
PYGIT=${PROJ}/scripts/

for i in tdu tpo;
do
 
    INPUT=$PROJ/empirical_bayesian_parents_input/ase_bayes_${i}_tdu_tpo_flag.csv
    OUTDIR=$PROJ/empirical_bayesian_parents_input/output/split_${i}_uo
    if [ ! -e $OUTDIR ]; then mkdir -p $OUTDIR; fi 

    # Split table
    python $PYGIT/splitTable.py -f $INPUT -o $OUTDIR --prefix split --header --nfiles 500
done

for i in tdu tpr;
do

    INPUT=$PROJ/empirical_bayesian_parents_input/ase_bayes_${i}_tdu_tpr_flag.csv
    OUTDIR=$PROJ/empirical_bayesian_parents_input/output/split_${i}_ur
    if [ ! -e $OUTDIR ]; then mkdir -p $OUTDIR; fi

    # Split table
    python $PYGIT/splitTable.py -f $INPUT -o $OUTDIR --prefix split --header --nfiles 500
done
