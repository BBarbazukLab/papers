##Job settings 
#PBS -N dros_norm 
#PBS -m abe 
##PBS -M ammorse@ufl.edu 
##Job Resources 
#PBS -l nodes=1:ppn=5
#PBS -l pmem=32gb 
#PBS -l walltime=160:00:00
#PBS -r n
#PBS -q bigmem

PROJ=/scratch/lfs/mcintyre/trago
OUT=/scratch/lfs/mcintyre/trago/dros_test_assembly/normalized_trinity
IN=/scratch/lfs/mcintyre/trago/dros_test_assembly/trimmomatic

module load trinity/r20140413
unset _JAVA_OPTIONS;
module load java/1.6.0_31
export _JAVA_OPTIONS="-Xms1g -Xmx30g"

## (1) normalize where 1st 10bp removed
normalize_by_kmer_coverage.pl \
        --seqType fq \
        --JM 190G \
        --max_cov 50 \
        --left_list $IN/r324_1_trim_paired.fq \
        --right_list $IN/r324_2_trim_paired.fq \
        --pairs_together \
        --PARALLEL_STATS \
        --JELLY_CPU 10 \
        --min_kmer_cov 2 \
        --output $OUT/r324_normalized_paired_reads

## (2) normalize where kept 10bp removed
normalize_by_kmer_coverage.pl \
        --seqType fq \
        --JM 190G \
        --max_cov 50 \
	--left_list $IN/r324_1_trim_with10_paired.fq \
        --right_list $IN/r324_2_trim_with10_paired.fq \
        --pairs_together \
        --PARALLEL_STATS \
        --JELLY_CPU 10 \
        --min_kmer_cov 2 \
        --output $OUT/r324_with10_normalized_paired_reads

