#!/bin/bash - 
#===============================================================================
# 
#   DESCRIPTION: Combine the split files from the ASE table 
# 
#===============================================================================

# Load python
module load python/2.7.6

# Set Variables
PROJ=/scratch/lfs/mcintyre/trago
PYGIT=/scratch/lfs/mcintyre/python.git

## only doing Tdu_tpr for Tms

for i in Tdu_Tdu_tpr Tpr_Tdu_tpr ;
do 

    INPUT=$PROJ/trago_output/emp_bayesian_parents/output/split_${i}

    OUTDIR=$PROJ/trago_output/emp_bayesian_parents/output
    if [ ! -e $OUTDIR ]; then mkdir -p $OUTDIR; fi 

    # combine table
    python $PYGIT/catTable.py -f $INPUT/split_*.csv --odir $OUTDIR --oname PG_emp_bayesian_results_${i}_parents.csv --header 

done
