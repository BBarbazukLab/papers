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
PROJ=/scratch/lfs/mcintyre/trago
PYGIT=/scratch/lfs/mcintyre/python.git

## only doing Tdu-tpr for Tms

for i in Tdu_Tdu_tpr Tpr_Tdu_tpr ;
do
 
    INPUT=$PROJ/trago_output/emp_bayesian_parents/input/ase_bayes_${i}_sbys_flag.csv

    OUTDIR=$PROJ/trago_output/emp_bayesian_parents/input/split_${i}
    if [ ! -e $OUTDIR ]; then mkdir -p $OUTDIR; fi 

    # Split table
    python $PYGIT/splitTable.py -f $INPUT -o $OUTDIR --prefix split --header --nfiles 500
done
