#!/usr/bin/env python

import sys, os, csv, argparse
csv.field_size_limit(1000000000)

##  This script modifies a tsv file
##  It reorders the columns in a tsv file that does not have fb# or chrom# removes all information after the initial pipe of the contig id, only leaving the 1st contig number


# Parse command line arguments
parser = argparse.ArgumentParser(description='clean up contig names in sam files')
parser.add_argument('-i','--input',dest='input_bed', action='store', required=True, help='input bed file [Required]')
parser.add_argument('-o','--out', dest='output_tsv', action='store', required=True, help='output tsv file [Required]')
args = parser.parse_args()


## open tsv file
with open(args.input_bed, 'rt') as bed_file:
    with open(args.output_tsv, "w") as output:
        for line in csv.reader(bed_file, delimiter='\t'):
            contig = line[0]
            start = line[1]
	    end = line[2]
            output.write(contig+"\t"+' '+"\t"+' '+"\t"+start+"\t"+end+"\n")
        
