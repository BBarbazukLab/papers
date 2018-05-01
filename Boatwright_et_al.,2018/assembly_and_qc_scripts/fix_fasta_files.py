#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import sys, os, re, csv

with open("/scratch/lfs/mcintyre/trago/references/test.fasta") as fa_file:
    with open("/scratch/lfs/mcintyre/trago/references/testout.fasta", "w") as output:
        for record in SeqIO.parse(fa_file, "fasta"):
            mult_ID = record.id
            ID = mult_ID.split("|")
            newID = ID[0]
            my_seqs =  SeqRecord(record.seq, id = newID, description='')

#            print my_seqs.id
            SeqIO.write(my_seqs, output, "fasta")
