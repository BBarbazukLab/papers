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
PYGIT=$PROJ/scripts

for i in Tm_tdu_tpo Tms_tdu_tpr;
do
 
    INPUT=$PROJ/empirical_bayesian_hybrids_input/ase_bayes_${i}_flag.csv
    OUTDIR=$PROJ/empirical_bayesian_hybrids_input/output/split_${i}
    if [ ! -e $OUTDIR ]; then mkdir -p $OUTDIR; fi 

    # Split table
    python $PYGIT/splitTable.py -f $INPUT -o $OUTDIR --prefix split --header --nfiles 500
done

