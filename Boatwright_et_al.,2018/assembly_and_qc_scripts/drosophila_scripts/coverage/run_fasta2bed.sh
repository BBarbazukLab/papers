#!/bin/bash - 
#===============================================================================
#
#          FILE: run_fb551_canonical_200bpJunction_fasta2bed.sh
# 
#         USAGE: ./run_fb551_canonical_200bpJunction_fasta2bed.sh 
# 
#   DESCRIPTION: An execution script to launch fasta2bed conversion
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: YOUR NAME (), 
#  ORGANIZATION: 
#       CREATED: 05/28/2013 12:54:32 PM EDT
#      REVISION:  ---
#===============================================================================

PROJ=$MCLAB/useful_dmel_data/flybase551
SCRIPTS=$PROJ/scripts
OUTPUT=$PROJ/output

FASTA=$OUTPUT/fb551_canonical_200bpJunctions.fasta
BED=$OUTPUT/fb551_canonical_200bpJunctions_fasta2bed.bed

python $SCRIPTS/fasta2bed.py \
    -f $FASTA \
    -o $BED 



