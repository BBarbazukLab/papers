#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import sys, os, re, itertools, operator, collections, csv, argparse
from itertools import islice

## script replaces spaces in fq read id with an underbar

# Parse command line arguments
#parser = argparse.ArgumentParser(description='clean up contig names in sam files')
#parser.add_argument('-i','--input',dest='input_fq', action='store', required=True, help='input fq file [Required]')
#parser.add_argument('-o','--out', dest='output_fq', action='store', required=True, help='fixed fq file [Required]')
#args = parser.parse_args()


with open('/scratch/lfs/mcintyre/trago/aln_test/list_R1.fastq', 'rb') as input_fq:
    with open('/scratch/lfs/mcintyre/trago/aln_test/testout.fastq', "w") as output_fq:
#        for line in input_fq:        
#            if input_fq.startswith('@'):
#                ID=re.sub(r' ', '_')
#                print ID
                    

        for record in SeqIO.parse(input_fq, "fastq-sanger"):
            space_ID = record.description
            ID = re.sub(r' ', '_', space_ID)
            my_seqs =  SeqRecord(record.seq, id = ID, description='', letter_annotations=record.letter_annotations)
#            print ID
            SeqIO.write(my_seqs, output_fq, "fastq")
