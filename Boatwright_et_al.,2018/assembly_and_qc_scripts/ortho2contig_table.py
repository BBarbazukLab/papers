#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import sys, os, re, csv, argparse

## script generates a table containing trago common orthos and their associated contig IDs


tdu_file = open('/scratch/lfs/mcintyre/trago/references/Tdu_common_ortho.fasta')
outTdu = open('/scratch/lfs/mcintyre/trago/references/orthoID2contigID_tdu.txt', "w")
for record in SeqIO.parse(tdu_file, "fasta"):
    mult_ID = record.id
    ID = mult_ID.split(",")
    orthoID = ID[0]
    tduID = ID[1]
    tdu_length = len(record)
    outTdu.write(orthoID+"\t"+tduID+"\t"+str(tdu_length)+"\n")

tpo_file = open('/scratch/lfs/mcintyre/trago/references/Tpo_common_ortho.fasta')
outTpo = open('/scratch/lfs/mcintyre/trago/references/orthoID2contigID_tpo.txt', "w")
for record in SeqIO.parse(tpo_file, "fasta"):
    mult_ID = record.id
    ID = mult_ID.split(",")
    orthoID = ID[0]
    tpoID = ID[1] 
    tpo_length = len(record)
    outTpo.write(orthoID+"\t"+tpoID+"\t"+str(tpo_length)+"\n")

tpr_file = open('/scratch/lfs/mcintyre/trago/references/Tpr_common_ortho.fasta')
outTpr = open('/scratch/lfs/mcintyre/trago/references/orthoID2contigID_tpr.txt', "w")         		    
for record in SeqIO.parse(tpr_file, "fasta"):
    mult_ID = record.id
    ID = mult_ID.split(",")
    orthoID = ID[0]
    tprID = ID[1]
    tpr_length = len(record)
    outTpr.write(orthoID+"\t"+tprID+"\t"+str(tpr_length)+"\n")
