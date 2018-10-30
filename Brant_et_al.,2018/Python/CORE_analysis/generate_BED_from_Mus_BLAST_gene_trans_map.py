# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 14:59:34 2016

@author: Lucas Boatwright
"""

from sys import argv, stderr
import argparse
from datetime import datetime

def parse_arguments():
    """Parse arguments passed to script"""
    parser = argparse.ArgumentParser(description=
    "This script is designed to generate a BED file using BLAST HSPs \
    \nand the filtered gene to transcript map.\n\n \
    \nUsage: python {0} -m reciprocated_blast_hits.gene_trans.txt \
-a parsed_Mus_to_Trinity.blast.txt\n".format(argv[0]),
    formatter_class = argparse.RawDescriptionHelpFormatter)

    requiredNamed = parser.add_argument_group('required arguments')

    requiredNamed.add_argument("-m", "--GENE_TRANS_MAP", type=str, required=True,\
    help="Reciprocated-best-hit gene transcript map.", action="store")
    
    requiredNamed.add_argument("-a", "--BLAST_A", type=str, required=True,\
    help="Parsed BLAST file from Mus to Acomys in tab-separated format.", 
    action="store")

    return parser.parse_args() 


def get_reciprocated_table(file_name):
    """Get reciprocated orthologs"""
    with open(file_name) as f:
        reciprocated_file = f.readlines()
    reciprocated_table = []
    for line in reciprocated_file:
        split_line = line.split()
        reciprocated_table.append((split_line[0], split_line[1]))
    return reciprocated_table    
        
        
def filter_blast_file(blast_file, reciprocated_table):
    """Filter blast file by reciprocal orthologs"""        
    filtered_blast = []
    with open(blast_file) as blast:
        for x, line in enumerate(blast):
            if ( (x+1)%10000 == 0 ):
                stderr.write("\tProcessed {0} lines.\n".format(str(x+1)))
            split_line = line.split()
            if ((split_line[0], split_line[2])) in reciprocated_table:
                filtered_blast.append(line)
                continue        
            else:
                stderr.write("\tWARNING: Unreciprocated BLAST hit: {0}-{1}".format(
                    split_line[0],split_line[2]))
    return filtered_blast
        
        
def write_filtered_blast(filtered_blast, blast_file_name):
    """Write filtered blast"""        
    output = open(blast_file_name + ".filtered", "w")
    for line in filtered_blast:
        output.write(line)
    output.close()
        

def generate_S_and_Q_BEDs(filtered_blast, blast_file_name):
    """Generate bed files for Queries and Subjects"""
    output = open(blast_file_name + ".filtered_S.bed", "w")
    # write A
    for line in filtered_blast:
        split_line = line.split()
        output.write("\t".join([split_line[2], 
            str(min(int(split_line[11])-1, int(split_line[12])-1)), 
            str(max(int(split_line[11])-1, int(split_line[12])-1)),
            ".",".", "+\n"]))
    output.close()
    
    output = open(blast_file_name + ".filtered_Q.bed", "w")
    #write B
    for line in filtered_blast:
        split_line = line.split()
        output.write("\t".join([split_line[0], 
            str(min(int(split_line[9])-1, int(split_line[10])-1)), 
            str(max(int(split_line[9])-1, int(split_line[10])-1)),
            ".",".", "+\n"]))
    output.close()
               
               
if __name__ == "__main__":
    start = datetime.now()
    args = parse_arguments()
    stderr.write("Executed: python {0} -m {1} -a {2}\n".format(argv[0],
                 args.GENE_TRANS_MAP,args.BLAST_A))    
    
    stderr.write("Reading gene transcript map.\n")
    reciprocated_table = get_reciprocated_table(args.GENE_TRANS_MAP)
    
    stderr.write("Reading parsed BLAST file.\n")   
    filtered_blast = filter_blast_file(args.BLAST_A, reciprocated_table)
    
    stderr.write("Writing filtered BLAST file.\n") 
    write_filtered_blast(filtered_blast, args.BLAST_A)
    
    stderr.write("Generating Subject and Query BED files.\n") 
    generate_S_and_Q_BEDs(filtered_blast, args.BLAST_A)
    
    stop = datetime.now()
    stderr.write("Runtime: {0}\n".format(str(stop - start)))
    
