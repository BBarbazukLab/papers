#!/bin/bash - 
#===============================================================================
# 
#   DESCRIPTION: Combine the split files from the ASE table 
# 
#===============================================================================

# Load python
module load python/2.7.6

# Set Variables
PROJ=/ufrc/barbazuk/lboat/Old_World_New_BED
PYGIT=/ufrc/barbazuk/lboat/Old_World_New_BED/scripts

for i in Tm_tdu_tpo Tms_tdu_tpr ;
do

    INPUT=$PROJ/empirical_bayesian_hybrids_output/split_${i}

    OUTDIR=$PROJ/empirical_bayesian_hybrids_output

    # combine table
    python $PYGIT/catTable.py -f $INPUT/split_*.csv --odir $OUTDIR --oname PG_emp_bayesian_results_${i}_hybrid.csv --header

done

