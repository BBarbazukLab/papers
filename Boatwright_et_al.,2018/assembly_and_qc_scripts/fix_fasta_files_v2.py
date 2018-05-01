#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import sys, os, re, csv, argparse

## script removes information from fasta id after the initial pipe, leaving only the 1st contig number

# Parse command line arguments
parser = argparse.ArgumentParser(description='clean up contig names in sam files')
parser.add_argument('-i','--input',dest='input_fa', action='store', required=True, help='input fa file [Required]')
parser.add_argument('-o','--out', dest='output_fa', action='store', required=True, help='cleaned up fa file [Required]')
args = parser.parse_args()


with open(args.input_fa) as fa_file:
    with open(args.output_fa, "w") as output:
        for record in SeqIO.parse(fa_file, "fasta"):
            mult_ID = record.id
            ID = mult_ID.split("|")
            newID = ID[0]
            my_seqs =  SeqRecord(record.seq, id = newID, description='')
            SeqIO.write(my_seqs, output, "fasta")
