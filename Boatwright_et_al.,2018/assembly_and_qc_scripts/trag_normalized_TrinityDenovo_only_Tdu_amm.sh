##Job settings
#PBS -N trin_on_norm
#PBS -m abe
##PBS -M ammorse@ufl.edu
##Job Resources
#PBS -l nodes=1:ppn=5
#PBS -l pmem=32gb
#PBS -l walltime=50:00:00
#PBS -r n
#PBS -q bigmem

PROJ=/scratch/lfs/mcintyre/trago
OUT=/scratch/lfs/mcintyre/trago/trago_output/trag_normalized_trinity_SC
IN=/scratch/lfs/mcintyre/trago/trago_output

module load trinity/r20130225
unset _JAVA_OPTIONS;
module load java/1.6.0_31
export _JAVA_OPTIONS="-Xms1g -Xmx30g"


## (1) Tdu normalized trinity

## rm -r $OUT/Tdu_dir_denovTrinityOutput
## mkdir $OUT/Tdu_dir_denovTrinityOutput

Trinity.pl \
	--seqType fq \
        --left $IN/left_Tdu.txt.normalized_K25_C50_pctSD100.fq \
        --right $IN/right_Tdu.txt.normalized_K25_C50_pctSD100.fq \
        --JM 90G \
        --bflyHeapSpaceInit 1G \
        --bflyHeapSpaceMax 30G \
        --CPU 5 \
        --output $OUT/Tdu_dir_denovTrinityOutput \
        --full_cleanup \
        --min_contig_length 200
