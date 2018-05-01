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

for i in tdu tpo ;
do

    INPUT=$PROJ/empirical_bayesian_parents_output/split_${i}_uo

    OUTDIR=$PROJ/empirical_bayesian_parents_output

    # combine table
    python $PYGIT/catTable.py -f $INPUT/split_*.csv --odir $OUTDIR --oname PG_emp_bayesian_results_${i}_parents_UO.csv --header

done

for i in tdu tpr ;
do

    INPUT=$PROJ/empirical_bayesian_parents_output/split_${i}_ur

    OUTDIR=$PROJ/empirical_bayesian_parents_output

    # combine table
    python $PYGIT/catTable.py -f $INPUT/split_*.csv --odir $OUTDIR --oname PG_emp_bayesian_results_${i}_parents_UR.csv --header

done
