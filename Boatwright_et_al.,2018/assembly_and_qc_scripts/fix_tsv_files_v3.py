#!/usr/bin/env python

import sys, os, csv, argparse
csv.field_size_limit(1000000000)

##  This script modifies a tsv file
##  It removes all information after the initial pipe of the contig id, only leaving the 1st contig number


# Parse command line arguments
parser = argparse.ArgumentParser(description='clean up contig names in sam files')
parser.add_argument('-i','--input',dest='input_tsv', action='store', required=True, help='input tsv file [Required]')
parser.add_argument('-o','--out', dest='output_tsv', action='store', required=True, help='cleaned up tsv file [Required]')
args = parser.parse_args()


## open sam file
with open(args.input_tsv, 'rt') as tsv_file:
    with open(args.output_tsv, "w") as output:
        for line in csv.reader(tsv_file, delimiter='\t'):
            mult_ID = line[0]
            ID=mult_ID.split("|")
            newID=ID[0]
            col1 = line[1]
	    col2 = line[2]	
	    col3 = line[3]
	    col4 = line[4]
            output.write(str(newID)+"\t"+col1+"\t"+col2+"\t"+col3+"\t"+col4+"\n")
        
