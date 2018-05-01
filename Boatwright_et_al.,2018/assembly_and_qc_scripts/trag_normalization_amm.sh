#!/bin/bash
#PBS -N Trin_norm
#PBS -M ammorse@ufl.edu
#PBS -m n
#PBS -r n
#PBS -o /scratch/lfs/mcintyre/trago/scripts/PBS_LOGS/
#PBS -j oe
#PBS -l nodes=1:ppn=10
#PBS -l walltime=200:00:00
#PBS -l pmem=24gb
#PBS -q bigmem


module load trinity/r20130225
unset _JAVA_OPTIONS;
module load java/1.6.0_31
export _JAVA_OPTIONS="-Xms1g -Xmx30g"

## Tdu normaliation WITHOUT 1st 10bp
normalize_by_kmer_coverage.pl \
	--seqType fq \
	--JM 190G \
	--max_cov 50 \
	--left_list /scratch/lfs/mcintyre/trago/trago_output/left_Tdu.txt \
	--right_list /scratch/lfs/mcintyre/trago/trago_output/right_Tdu.txt \
	--pairs_together \
	--PARALLEL_STATS \
	--JELLY_CPU 10 \
	--min_kmer_cov 2 \
	--output /scratch/lfs/mcintyre/trago/trago_output/trago_norm/Tdu_normalized_paired_reads

## Tdu normaliation WITH 1st 10bp 
normalize_by_kmer_coverage.pl \
        --seqType fq \
        --JM 190G \
        --max_cov 50 \
        --left_list /scratch/lfs/mcintyre/trago/trago_output/left_Tdu_10.txt \
        --right_list /scratch/lfs/mcintyre/trago/trago_output/right_Tdu_10.txt \
        --pairs_together \
        --PARALLEL_STATS \
        --JELLY_CPU 10 \
        --min_kmer_cov 2 \
        --output /scratch/lfs/mcintyre/trago/trago_output/trago_norm/Tdu_10_normalized_paired_reads


## Tpo normaliation WITHOUT 1st 10bp
normalize_by_kmer_coverage.pl \
        --seqType fq \
        --JM 190G \
        --max_cov 50 \
        --left_list /scratch/lfs/mcintyre/trago/trago_output/left_Tpo.txt \
        --right_list /scratch/lfs/mcintyre/trago/trago_output/right_Tpo.txt \
        --pairs_together \
        --PARALLEL_STATS \
        --JELLY_CPU 10 \
        --min_kmer_cov 2 \
        --output /scratch/lfs/mcintyre/trago/trago_output/trago_norm/Tpo_normalized_paired_reads

## Tpo normaliation WITH 1st 10bp 
normalize_by_kmer_coverage.pl \
        --seqType fq \
        --JM 190G \
        --max_cov 50 \
        --left_list /scratch/lfs/mcintyre/trago/trago_output/left_Tpo_10.txt \
        --right_list /scratch/lfs/mcintyre/trago/trago_output/right_Tpo_10.txt \
        --pairs_together \
        --PARALLEL_STATS \
        --JELLY_CPU 10 \
        --min_kmer_cov 2 \
        --output /scratch/lfs/mcintyre/trago/trago_output/trago_norm/Tpo_10_normalized_paired_reads

## Tpr normaliation WITHOUT 1st 10bp
normalize_by_kmer_coverage.pl \
        --seqType fq \
        --JM 190G \
        --max_cov 50 \
        --left_list /scratch/lfs/mcintyre/trago/trago_output/left_Tpr.txt \
        --right_list /scratch/lfs/mcintyre/trago/trago_output/right_Tpr.txt \
        --pairs_together \
        --PARALLEL_STATS \
        --JELLY_CPU 10 \
        --min_kmer_cov 2 \
        --output /scratch/lfs/mcintyre/trago/trago_output/trago_norm/Tpr_normalized_paired_reads

## Tpr normaliation WITH 1st 10bp
normalize_by_kmer_coverage.pl \
        --seqType fq \
        --JM 190G \
        --max_cov 50 \
        --left_list /scratch/lfs/mcintyre/trago/trago_output/left_Tpr_10.txt \
        --right_list /scratch/lfs/mcintyre/trago/trago_output/right_Tpr_10.txt \
        --pairs_together \
        --PARALLEL_STATS \
        --JELLY_CPU 10 \
        --min_kmer_cov 2 \
        --output /scratch/lfs/mcintyre/trago/trago_output/trago_norm/Tpr_10_normalized_paired_reads
