##Job settings 
#PBS -N Tdu_sub.sh 
#PBS -o Tdu_sub.out 
#PBS -e Tdu_sub.err 
#PBS -m abe 
##PBS -M schamala@ufl.edu 
##Job Resources 
#PBS -l nodes=1:ppn=5
#PBS -l pmem=32gb 
#PBS -l walltime=160:00:00
#PBS -r n
#PBS -q bigmem

#PBS -o /scratch/lfs/mcintyre/trago/Tragopogon_Normalized_Trinity_SC 


module load trinity/r20130225
unset _JAVA_OPTIONS;
module load java/1.6.0_31
export _JAVA_OPTIONS="-Xms1g -Xmx30g"

rm -r Tdu_dir_denovTrinityOutput
mkdir Tdu_dir_denovTrinityOutput

/apps/trinity/r20130225/util/../Trinity.pl --seqType fq --left Tdu_left_reads_list.txt.normalized_K25_C50_pctSD100.fq --right Tdu_right_reads_list.txt.normalized_K25_C50_pctSD100.fq --JM 90G --bflyHeapSpaceInit 1G --bflyHeapSpaceMax 30G --CPU 5 --output Tdu_dir_denovTrinityOutput --full_cleanup --min_contig_length 200  

rm -r Tpr_dir_denovTrinityOutput
mkdir Tpr_dir_denovTrinityOutput

/apps/trinity/r20130225/util/../Trinity.pl --seqType fq --left Tpr_left_reads_list.txt.normalized_K25_C50_pctSD100.fq --right Tpr_right_reads_list.txt.normalized_K25_C50_pctSD100.fq --JM 90G --bflyHeapSpaceInit 1G --bflyHeapSpaceMax 30G --CPU 5 --output Tpr_dir_denovTrinityOutput --full_cleanup --min_contig_length 200  


rm -r Tpo_dir_denovTrinityOutput
mkdir Tpo_dir_denovTrinityOutput

/apps/trinity/r20130225/util/../Trinity.pl --seqType fq --left Tpo_left_reads_list.txt.normalized_K25_C50_pctSD100.fq --right Tpo_right_reads_list.txt.normalized_K25_C50_pctSD100.fq --JM 90G --bflyHeapSpaceInit 1G --bflyHeapSpaceMax 30G --CPU 5 --output Tpo_dir_denovTrinityOutput --full_cleanup --min_contig_length 200  


