#!/usr/bin/env python

import sys, os, csv, argparse
csv.field_size_limit(1000000000)

# Parse command line arguments
parser = argparse.ArgumentParser(description='clean up contig names in sam files')
parser.add_argument('-i','--input',dest='input_sam', action='store', required=True, help='input sam file [Required]')
parser.add_argument('-o','--out', dest='output_sam', action='store', required=True, help='cleaned up sam file [Required]')
args = parser.parse_args()


## open sam file
with open(args.input_sam, 'rt') as sam_file:
    with open(args.output_sam, "w") as output:
        for line in csv.reader(sam_file, delimiter='\t'):
            col0 = line[0]
            col1 = line[1]
            mult_ID=line[2]
            col3=line[3]
            col4=line[4]
            col5=line[5]
            col6=line[6]
            col7=line[7]
            col8=line[8]
            col9=line[9]
            col10=line[10]
            col11=line[11]
            col12=line[12]
            col13=line[13]
	    mult_ID=line[2]
	    ID=mult_ID.split("|")
	    newID=ID[0]
            output.write(col0 + "\t" + col1 + "\t" + str(newID) + "\t" + col3 + "\t" +col4 + "\t" +col5 + "\t" +col6 + "\t" +col7 + "\t" +col8 + "\t" +col9 + "\t" +col10 + "\t" +col11 + "\t" +col12 + "\t" +col13+"\n")
        
