

#Ran on our local machine

/apps/trinity/r20130225/util/normalize_by_kmer_coverage.pl --seqType fq --JM 190G --max_cov 50 --left_list Tdu_left_reads_list.txt --right_list Tdu_right_reads_list.txt --pairs_together --PARALLEL_STATS --JELLY_CPU 16 --min_kmer_cov 2 --output Tdu_normalized_paired_reads

/apps/trinity/r20130225/util/normalize_by_kmer_coverage.pl --seqType fq --JM 190G --max_cov 50 --left_list Tpr_left_reads_list.txt --right_list Tpr_right_reads_list.txt --pairs_together --PARALLEL_STATS --JELLY_CPU 16 --min_kmer_cov 2 --output Tpr_normalized_paired_reads 


/apps/trinity/r20130225/util/normalize_by_kmer_coverage.pl --seqType fq --JM 190G --max_cov 50 --left_list Tpo_left_reads_list.txt --right_list Tpo_right_reads_list.txt --pairs_together --PARALLEL_STATS --JELLY_CPU 16 --min_kmer_cov 2 --output Tpo_normalized_paired_reads 

