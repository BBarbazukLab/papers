#!/usr/bin/env python

import logging
import argparse
from Bio import SeqIO
import os.path

def getOptions():
    """ Function to pull in arguments """
    parser = argparse.ArgumentParser(description="Reads a FASTQ file and outputs sequences greater than 1000bp.")
    parser.add_argument("-i", "--input", dest="fname", action='store', required=True, help="Name of input FASTQ file [Required]")
    parser.add_argument("-o", "--out", dest="out", action='store', required=True, help="Name of output csv file [Required]")
    parser.add_argument("-l", "--length", dest="length", action='store', required=True, help="truncate at this length [Required]")
    args = parser.parse_args()
    return(args)


def fastq_reader(args):
    """ Function to read FASTQ file and parse the sequences """
    long_seqs = []
    for record in SeqIO.parse(args.fname, "fasta") :
        if len(record.seq) > args.length :
            long_seqs.append(record)

    print "Found %i long sequences in fname" % len(long_seqs)

    output_handle = open(args.out, "w")
    SeqIO.write(long_seqs, output_handle, "fasta")
    output_handle.close()


def main():
    """ MAIN Function to execute everything """ 
    args = getOptions()
    fastq_reader(args)

if __name__=='__main__':
    main()

