##Job settings 
#PBS -N trin_dros 
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


## (1) run where 1st 10bp removed
##Trinity --normalize_reads --seqType fq \
##--left $IN/r324_1_trim_paired.fq \
##--right $IN/r324_2_trim_paired.fq \
##--JM 90G \
##--bflyHeapSpaceInit 1G \
##--bflyHeapSpaceMax 30G \
##--CPU 5 \
##--output $OUT/r324_paired_trinityOutput \
##--full_cleanup \
##--min_contig_length 200  


## (2) run where left 1st 10bp on 
Trinity --normalize_reads --seqType fq \
--left $IN/r324_1_trim_with10_paired.fq \
--right $IN/r324_2_trim_with10_paired.fq \
--JM 90G \
--bflyHeapSpaceInit 1G \
--bflyHeapSpaceMax 30G \
--CPU 5 \
--output $OUT/r324_paired_with10_trinityOutput \
--full_cleanup \
--min_contig_length 200

