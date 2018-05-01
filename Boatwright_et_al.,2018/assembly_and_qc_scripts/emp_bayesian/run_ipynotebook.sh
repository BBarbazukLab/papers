#!/bin/bash
PROJ=$MCLAB/cegs_ase_paper
DOCS=$PROJ/documentation
SCRIPTS=$PROJ/scripts/emp_bayesian

# Run the IPython Notebook
cd $SCRIPTS
ipython notebook
