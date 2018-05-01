# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 11:16:34 2016

@author: Lucas Boatwright
"""

from sys import argv, stderr, exit
import argparse
from datetime import datetime

def parse_arguments():
    """Parse arguments passed to script"""
    parser = argparse.ArgumentParser(description=
    "This script is designed to identify reciprocal best-hits given two \
    \nparsed BLAST files. Output includes reciprocated_blast_hits.txt and \
    \ntwo BED files representing HSP coordinates from BLAST_A for query \
    \nand subject.\n\
    \n\tUsage: python {0} -a A_to_B_BLAST -b B_to_A_BLAST\
    \n".format(
    argv[0]), formatter_class = argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-a", "--BLAST_A", type=str, required=True,
    help="Parsed BLAST file from species A in tab-separated format.", 
        action="store")
        
    parser.add_argument("-b", "--BLAST_B", type=str, required=True,
    help="Parsed BLAST file from species B in tab-separated format.", 
        action="store")
        
    parser.add_argument("-p", "--PROGRAM", type=str, required=True,
    help="Program used to generate tabular blast -- NCBI or WU parsed.", 
        action="store")

    return parser.parse_args() 


def get_dict_from_blast(file_name, subject_index):
    """Generate a dictionary from tabular BLAST (one-one) output"""
    blast_dict = {}
    with open(file_name) as file:
        for line in file:
            split_line = line.split()
            blast_dict[split_line[0]] = split_line[subject_index]
    return blast_dict

    
def check_reciprocity(dict_a, dict_b):
    """Check to see BLAST homology is reciprocated""" 
    reciprocated = {}
    for key, value in dict_a.items():               
        try:
            if dict_b[value]==key:
                reciprocated[key] = value
            else:
                continue
        except KeyError:
            continue
    output = open("reciprocated_blast_hits.txt","w")
    for key, value in reciprocated.items():
        if value != []:
            output.write(key + "\t" + value + "\n")
            output.flush()
    output.close()


def get_reciprocated_table():
    """Get reciprocated orthologs"""
    with open("reciprocated_blast_hits.txt") as f:
        reciprocated_file = f.readlines()
    reciprocated_table = []
    for line in reciprocated_file:
        split_line = line.split()
        reciprocated_table.append((split_line[0], split_line[1]))
    return reciprocated_table    
        
        
def filter_blast_file(blast_file, reciprocated_table, subject_index):
    """Filter blast file by reciprocal orthologs"""       
    filtered_blast = []
    with open(blast_file) as blast:
        for x, line in enumerate(blast):
            if ( (x+1)%10000 == 0 ):
                stderr.write("\tProcessed {0} lines.\n".format(str(x+1)))
            split_line = line.split()
            if ((split_line[0], split_line[subject_index])) in reciprocated_table:
                filtered_blast.append(line)
    return filtered_blast
        

def generate_S_and_Q_BEDs(filtered_blast, blast_file_name, subject_index):
    """Generate bed files for Queries and Subjects"""
    if subject_index == 1:
        q_start = 6
        q_stop = 7
        s_start = 8
        s_stop = 9
    elif subject_index == 2:
        q_start = 9
        q_stop = 10
        s_start = 11
        s_stop = 12
    output = open(blast_file_name + ".filtered_S.bed", "w")
    # write B
    for line in filtered_blast:
        split_line = line.split()
        output.write("\t".join([split_line[subject_index], 
            str(min(int(split_line[s_start])-1, int(split_line[s_stop])-1)), 
            str(max(int(split_line[s_start])-1, int(split_line[s_stop])-1)),
            ",".join([split_line[0],split_line[subject_index]]),".", "+\n"]))
    output.close()
    
    output = open(blast_file_name + ".filtered_Q.bed", "w")
    #write A
    for line in filtered_blast:
        split_line = line.split()
        output.write("\t".join([split_line[0], 
            str(min(int(split_line[q_start])-1, int(split_line[q_stop])-1)), 
            str(max(int(split_line[q_start])-1, int(split_line[q_stop])-1)),
            ",".join([split_line[0],split_line[subject_index]]),".", "+\n"]))
    output.close()
    
    
if __name__ == "__main__":
    start = datetime.now()
    args = parse_arguments()
    if args.PROGRAM.lower() == "ncbi":
        subject_index = 1
    elif args.PROGRAM.lower() == "wu":
        subject_index = 2
    else:
        stderr.write("Invalid program option. Use NCBI or WU.\n")
        exit(1)
    stderr.write("Generating dictionary from Species A BLAST to B.\n")
    dict_a = get_dict_from_blast(args.BLAST_A, subject_index)
    stderr.write("Generating dictionary from Species B BLAST to A.\n")
    dict_b = get_dict_from_blast(args.BLAST_B, subject_index)
    stderr.write("Identifying reciprocated hits.\n")
    check_reciprocity(dict_a, dict_b)
    stderr.write("Generating BED files from gene-transcript map.\n")
    reciprocated_table = get_reciprocated_table()
    filtered_blast = filter_blast_file(args.BLAST_A, reciprocated_table, 
                                       subject_index)
    generate_S_and_Q_BEDs(filtered_blast, args.BLAST_A, subject_index)
    stop = datetime.now()
    stderr.write("Runtime: {0}".format(str(stop - start)))