#!/bin/bash

module add bedtools
module add python

bedtools getfasta -bed TDU-tpo_overlaps_WRT_orthologs.bed -fi TDU_consensed_contigs_500-15k.fasta -fo bed_filtered_TDU-tpo.fasta
bedtools getfasta -bed TDU-tpr_overlaps_WRT_orthologs.bed -fi TDU_consensed_contigs_500-15k.fasta -fo bed_filtered_TDU-tpr.fasta

bedtools getfasta -bed TPO-tdu_overlaps_WRT_orthologs.bed -fi TPO_consensed_contigs_500-15k.fasta -fo bed_filtered_TPO-tdu.fasta
bedtools getfasta -bed TPR-tdu_overlaps_WRT_orthologs.bed -fi TPR_consensed_contigs_500-15k.fasta -fo bed_filtered_TPR-tdu.fasta

pyfasta info --gc -n -1 bed_filtered_TDU-tpo.fasta > bed_filtered_TDU-tpo.stats
pyfasta info --gc -n -1 bed_filtered_TDU-tpr.fasta > bed_filtered_TDU-tpr.stats

pyfasta info --gc -n -1 bed_filtered_TPO-tdu.fasta > bed_filtered_TPO-tdu.stats
pyfasta info --gc -n -1 bed_filtered_TPR-tdu.fasta > bed_filtered_TPR-tdu.stats

grep ">" bed_filtered_TDU-tpo.stats | cut -d ':' -f 1,4 | sed 's/:/,/g' | sed 's/%//g' > bed_filtered_TDU-tpo.stats.csv
grep ">" bed_filtered_TDU-tpr.stats | cut -d ':' -f 1,4 | sed 's/:/,/g' | sed 's/%//g' > bed_filtered_TDU-tpr.stats.csv

grep ">" bed_filtered_TPO-tdu.stats | cut -d ':' -f 1,4 | sed 's/:/,/g' | sed 's/%//g' > bed_filtered_TPO-tdu.stats.csv
grep ">" bed_filtered_TPR-tdu.stats | cut -d ':' -f 1,4 | sed 's/:/,/g' | sed 's/%//g' > bed_filtered_TPR-tdu.stats.csv

