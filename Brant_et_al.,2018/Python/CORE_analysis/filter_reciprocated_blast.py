# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 12:11:26 2016

@author: Lucas Boatwright
"""

from sys import argv, stderr
import argparse
from datetime import datetime

def parse_arguments():
    """Parse arguments passed to script"""
    parser = argparse.ArgumentParser(description=
    "This script is designed to filter reciprocated best hits based upon a \
    \nprovided quality threshold.\n\n \
    \nUsage: python {0} -m reciprocated_blast_hits.gene_trans.txt \
-a parsed_Mus_to_Trinity.blast.txt\n".format(argv[0]),
    formatter_class = argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-m", "--GENE_TRANS_MAP", type=str, required=True,\
    help="Reciprocated-best-hit gene transcript map", action="store")
    
    parser.add_argument("-a", "--BLAST_A", type=str, required=True,\
    help="Parsed BLAST file from Mus to Acomys in tab-separated format.", 
    action="store")

    return parser.parse_args() 
        

def open_blast(file_name):
    with open(file_name) as f:
        parsed_blast = f.readlines()
    return parsed_blast


def hits_to_dict(file_name):
    recip_hits = {}
    with open(file_name) as f:
        for line in f:
            list = line.split()
            recip_hits[(list[0].rstrip(), list[1].rstrip())]=1
    return recip_hits


def search_reciprocated_hits(recip_hits, parsed_blast, blast_file_name):
    output = open(blast_file_name.split(".txt")[0] + ".filtered.txt","w")
    for line in parsed_blast:
        split_line = line.split()
        try:
            if (recip_hits[(split_line[0].rstrip(), split_line[2].rstrip())]==1):
                if float(split_line[5]) <= 1e-5:
                    output.write(line)
                    output.flush()
        except KeyError:
            stderr.write("\t" + split_line[0] + " and " + split_line[2] + 
                " not reciprocated\n")
    output.close()

    
if __name__ == "__main__":
    start = datetime.now()
    args = parse_arguments()
    stderr.write("Executed: python {0} -m {1} -a {2}\n".format(argv[0],
                 args.GENE_TRANS_MAP,args.BLAST_A))    
    
    stderr.write("Reading parsed BLAST file.\n")
    parsed_blast = open_blast(args.BLAST_A)
    
    stderr.write("Reading gene transcript map.\n")
    recip_hits = hits_to_dict(args.GENE_TRANS_MAP)
    
    stderr.write("Filtering BLAST file.\n")
    search_reciprocated_hits(recip_hits, parsed_blast, args.BLAST_A)
    
    stop = datetime.now()
    stderr.write("Runtime: {0}\n".format(str(stop - start)))

