# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 12:28:24 2017

@author: Lucas Boatwright
"""

from sys import argv, stderr
import argparse
from datetime import datetime

def parse_arguments():
    """Parse arguments passed to script"""
    parser = argparse.ArgumentParser(description=
    "This script is designed to take a parsed BLAST file and generate a \
    \ntable of groups. The grouping method is greedy.\n\n \
    Example: python {0} -b parsed_BLAST".format(argv[0]),
    formatter_class = argparse.RawDescriptionHelpFormatter)
    
    requiredNamed = parser.add_argument_group('required arguments')

    requiredNamed.add_argument("-b", "--BLAST", type=str, required=True,\
    help="Parsed BLAST file used for grouping transcripts.", action="store")

    return parser.parse_args() 
        
        
def group_self_blast(file_name):
    """Generate groups for parsed self-BLAST files."""
    # Initialize lists to be used
    temp_array = []     # current group
    group_array = []    # list of groups
    num = 1             # line number -- only used for first
    all_contigs=[]      # all contigs that have been checked
    
    with open(file_name, 'r') as file:
        for line in file:
            split_line = line.rsplit()
            if num==1:
                # Start group at first split_line -- not used again
                current_query = split_line[0]
                temp_array.append(split_line[0])
                temp_array.append(split_line[2])
                all_contigs.append(split_line[0])
                all_contigs.append(split_line[2])
                num += 1
            elif ((split_line[0]==current_query) and (split_line[2] not in all_contigs)):
                # If still on same query and subject is new -- add it to group
                temp_array.append(split_line[2])
                all_contigs.append(split_line[2])
            elif ((split_line[0] in all_contigs) and (split_line[2] in all_contigs)):
                # If both query and subject are accounted for -- continue
                continue
            elif ((split_line[0] in all_contigs) and (split_line[2] not in all_contigs)):
                # If new subject, add it to the appropriate array
                for i in group_array:
                    if split_line[0] in i:
                        i.append(split_line[2])
                all_contigs.append(split_line[2])
            elif ((split_line[0] not in all_contigs) and (split_line[2] in all_contigs)):
                # If new query but subject is grouped, add to appropriate array
                for i in group_array:
                    if split_line[2] in i:
                        i.append(split_line[0])
                all_contigs.append(split_line[0])
            else:
                # If new query -- add temp_array to group array
                # and start process over for new query
                group_array.append(temp_array)
                temp_array = []
                # set current query
                current_query = split_line[0]
                temp_array.append(split_line[0])
                all_contigs.append(split_line[0])
                
                # if transcript is new, add to list
                if(split_line[2] not in all_contigs):
                    temp_array.append(split_line[2])
                    all_contigs.append(split_line[2])

    # write groups to file
    output = open("grouped_hits.txt",'w')
    for i in group_array:
        output.write(" ".join(i) + "\n")
    output.close()
    
    
if __name__ == "__main__":
    start = datetime.now()
    args = parse_arguments()
    stderr.write("Executed: python {0} -b {1}\n".format(argv[0],args.BLAST))
    
    stderr.write("Grouping transcripts.\n")
    group_self_blast(args.BLAST)    
    
    stop = datetime.now()
    stderr.write("Runtime: {0}\n".format(str(stop - start)))
