#!/bin/bash - 
#===============================================================================
#
#         USAGE: ./sync_hpc_xsub_3ci_dna_to_rlts.sh 
# 
#   DESCRIPTION: This script syncs the contents of a folder on HPC scratch to this storage location.
#                order from this folder to that folder   
#===============================================================================

rsync -rulvz --exclude './sync_hpc_xsub_3ci_dna_to_rlts.sh' /scratch/lfs/mcintyre/Xsub_3ci_genome_project/original_data/DNA/  /rlts/mcintyre/Xsub_3ci_project_DNA/
