#!/usr/bin/env python
import argparse
import itertools
import collections
from Bio.SeqIO.QualityIO import FastqGeneralIterator as FGI


# Parse command line arguments
parser = argparse.ArgumentParser(description='Counting total and unique fqs')
parser.add_argument('-i','--input',dest='fq', action='store', required=True, help='A list of fq file [Required]')
parser.add_argument('-o','--out', dest='out', action='store', required=True, help='Output file for counts in csv format [Required]')
args = parser.parse_args()


def dict():
    with open(args.fq,'r') as FQ:
        seqs = itertools.islice(FQ,1,None,4)
        counts = collections.Counter(seqs)
        return counts

counts = dict()

def counting():
    with open(args.out,'w') as dataout:
        with open(args.fq,'r') as FQ:
            for (header,seq,qual) in FGI(FQ):
                currCount = counts[seq + "\n"]
                short_head = header.split(' ')[0]
                if currCount > 1:
                    flag = short_head + "," + "1"
                else:
                    flag = short_head + "," + "0"
                myout = [flag]
                dataout.write(','.join(str(x) for x in myout) + "\n")

counting()
