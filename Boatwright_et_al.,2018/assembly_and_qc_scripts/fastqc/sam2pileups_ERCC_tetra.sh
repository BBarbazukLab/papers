#!/bin/bash
#PBS -M polvadore@ufl.edu
#PBS -m n
#PBS -q bio
#PBS -r n
#PBS -j oe
#PBS-o /bio/mcintyre/trago/scripts/fastqc/PBS_LOGS
#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=4
#PBS -l pmem=5gb
#PBS -t 1-36

module load samtools/0.1.18

### set directories
PROJ=/bio/mcintyre/trago
ALNS=$PROJ/ercc_alignments/tetraploid
REF=/bio/mcintyre/references/ERCC_Ambion

#### I am using an SGE Array (pulls out each name from the file and calls it 'design')
DESIGN_FILE=$PROJ/outfiles/fastqc/tetraploids/design_files/all_trago_tetra.txt
DESIGN=$(cat $DESIGN_FILE | head -n $PBS_ARRAYID | tail -n 1)  #steps through file_list.txt by row
ID=$(basename "$DESIGN")

### convert SAM to BAM
samtools view -ut $REF/ERCC92.fa.fai -o $ALNS/aln_ERCC_${ID}.bam $ALNS/aln_ERCC_${ID}.sam 

### sort BAM
samtools sort $ALNS/aln_ERCC_${ID}.bam $ALNS/aln_sorted_ERCC_${ID}

### create pileup file
samtools mpileup -f $REF/ERCC92.fa.fai $ALNS/aln_sorted_ERCC_${ID}.bam > $ALNS/pileups/${ID}_ERCC.pileup
