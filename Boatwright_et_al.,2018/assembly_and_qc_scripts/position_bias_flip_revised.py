#!/usr/bin/python

import logging
import argparse
from Bio import SeqIO
import os.path

def getOptions():
    """ Function to pull in arguments """
    parser = argparse.ArgumentParser(description="Reads a FASTQ file and calculates position bias.")
    parser.add_argument("-i", "--input", dest="fname", action='store', required=True, help="Name of input FASTQ file [Required]")
    parser.add_argument("-o", "--out", dest="out", action='store', required=True, help="Name of output csv file [Required]")
    parser.add_argument("-g", "--log", dest="log", action='store', required=False, help="Log File")
    args = parser.parse_args()
    return(args)

def setLogger(fname,loglevel):
    """ Function to handle error logging """
    logging.basicConfig(filename=fname, level=loglevel, format='%(asctime)s - %(levelname)s - %(message)s')




def fastq_reader(fname):
    """ Function to read fastq file and parse the sequences """
    with open(fname, "rU") as handle:
        file=list()
        for record in SeqIO.parse(handle, "fastq"):
            file.append(record.seq)
    return(file)

def position(x):  
    """ Function to calculate position bias """           
    position=list()
    num_sequences = len(x)    ## this is number of different sequences 
    flip = zip(*x)  ## flip data 90 degrees
    num_readbases = len(flip)    ## this is the number of readbases for each sequence
    for index, item in enumerate(flip):
        G_count = float(item.count('G')) # count the number of Gs, etc
        G = (G_count / num_sequences)    # convert to fraction of total
        A_count = float(item.count('A'))
        A = (A_count / num_sequences)
        T_count = float(item.count('T'))
        T = (T_count / num_sequences)
        C_count = float(item.count('C')) 
        C = (C_count / num_sequences)
    	results=(index, G, A, T, C)
     	position.append(results)
    return(position) 

def output_csv(outname,y):
    """ Function to output results to CSV """
    header = 'position G A T C'.split()
    with open(outname,'w') as f:
        f.write(','.join(str(x) for x in header) + '\n')
        for tup in y:
            line = ','.join(str(x) for x in tup)
            f.write(line + '\n')


def main(): 
    """ MAIN Function to execute everything """
    # Turn on Logging if option -g was given
    args = getOptions()\
    
    if args.log:
        setLogger(args.log,logging.INFO)

    logging.info("Converting '%s' to '%s'" % (args.fname,args.out))
    out = os.path.abspath(args.out)

    file_in = fastq_reader(args.fname)
    file_out = position(file_in)
    output_csv(args.out,file_out)
    logging.info("Finished converting '%s' to '%s'" % (args.fname,args.out))


if __name__=='__main__':
    main()
    logging.info("Script complete.")            

