# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 11:16:34 2016

@author: Lucas Boatwright
"""

from sys import argv, stderr
import argparse
from datetime import datetime

def parse_arguments():
    """Parse arguments passed to script"""
    parser = argparse.ArgumentParser(description=
    "This script is designed to identify reciprocal-best-hits given two \
    \nparsed BLAST files. Output includes reciprocated_blast_hits.txt and \
    \nreciprocated_blast_hits.gene_trans.txt\n\
    \n\tUsage: python {0} -a mus_to_aco_parsed_BLAST -b aco_to_mus_parsed_BLAST\
    \n\n\
    \nExample BLAST_A (sequence may hit multiple targets): \n\
    \nrna10340	101674	Trinity_contig_252389	66445	282021	0.	91.9	91.9901697537887	65919	98486	32576	58	65975	-\
    \nrna10340	101674	Trinity_contig_252388	66445	1621	0.	80.4	80.4878048780488	533	32644	32124	65919	66445	-\
    \nrna111265	843	Trinity_contig_22320	1097	1955	5.4e-106	83.1	83.1010452961672	574	582	14	520	1091	-\
    \n\nExample BLAST_B (one sequence to one sequence): \
    \nTrinity_contig_252388	66865	rna10340	101674	280762	0.	91.6	91.6699410609037	66170	66192	58	32362	98486	-\
    \nTrinity_contig_158768	524	rna90947	978	276	4.6e-05	79.5	79.5698924731183	93	339	429	289	380	+\
    \nTrinity_contig_167508	524	rna45620	4575	179	0.99993	61	61.0655737704918	244	289	516	1129	1357	+".format(
    argv[0]), formatter_class = argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-a", "--BLAST_A", type=str, required=True,\
    help="Parsed BLAST file from species A in tab-separated format.", 
        action="store")
        
    parser.add_argument("-b", "--BLAST_B", type=str, required=True,\
    help="Parsed BLAST file from species B in tab-separated format.", 
        action="store")

    return parser.parse_args() 


def get_dictionary_from_mus(file_name):
    """Generate a dictionary from the mus blast (one-many)"""
    mus_dict = {}
    with open(file_name) as f:
        for line in f:
            split_line = line.split()
            try:
                temp_values = mus_dict.pop(split_line[0])
                temp_values.append(split_line[2])
                mus_dict[split_line[0]] = temp_values
            except KeyError:
                mus_dict[split_line[0]] = [split_line[2]]
    for key, value in mus_dict.items():
        if len(value) > 1:
            mus_dict[key] = list(set(value))
    return mus_dict


def get_dictionary_from_acomys(file_name):
    """Generate a dictionary from the acomys blast (one-one)"""
    aco_dict = {}
    with open(file_name) as aco:
        for line in aco:
            split_line = line.split()
            aco_dict[split_line[0]] = split_line[2]
    return aco_dict

    
def check_reciprocity(mus_dict, aco_dict):
    """Check to see trinity-mus homology is reciprocated""" 
    for key, values in mus_dict.items():
        new_values = []
        for i in values:
            try:
                if aco_dict[i]==key:
                    new_values.append(i)
                else:
                    continue
            except KeyError:
                continue
        mus_dict[key] = new_values
    output = open("reciprocated_blast_hits.txt","w")
    for key, value in mus_dict.items():
        if value != []:
            output.write(key + "\t" + ",".join(value) + "\n")
            output.flush()
    output.close()
    
    
def get_reciprocated_gene_trans_map():
    """Generate a gene to transcript map from the reciprocated BLAST file"""
    output = open("reciprocated_blast_hits.gene_trans.txt",'w')
    with open("reciprocated_blast_hits.txt") as file:
        for line in file:
            split_line = line.split()
            for tr in split_line[1].split(","):
                output.write(split_line[0] + "\t" + tr + "\n")
                output.flush()
    output.close()
    
    
if __name__ == "__main__":
    start = datetime.now()
    args = parse_arguments()
    stderr.write("Generating dictionary from Mus BLAST to Acomys.\n")
    mus_dict = get_dictionary_from_mus(args.BLAST_A)
    stderr.write("Generating dictionary from Acomys BLAST to Mus.\n")
    aco_dict = get_dictionary_from_acomys(args.BLAST_B)
    stderr.write("Identifying reciprocated hits.\n")
    check_reciprocity(mus_dict, aco_dict)
    stderr.write("Generating reciprocated gene to transcript map.\n")
    get_reciprocated_gene_trans_map()
    stop = datetime.now()
    stderr.write("Runtime: {0}".format(str(stop - start)))

    
    
