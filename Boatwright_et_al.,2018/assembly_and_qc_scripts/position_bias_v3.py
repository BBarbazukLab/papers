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
    """ Function to read FASTQ file and parse the sequences """
    with open(fname, "rU") as handle:
        file=list()
        for record in SeqIO.parse(handle, "fastq"):
            file.append(record.seq)
    return(file)

def position(x):  
    """ Function to calculate position bias """
    position=list()
    num_seqs = len(x)
    num_bases = len(x[1])
    for j in range(0,num_bases): #each base
        count_A=0 ; count_T=0 ; count_C=0 ; count_G=0 ; count_other=0 ; total=0
        for i in range(0,num_seqs):  #each sequence
            if x[i][j]=='A':
                count_A=count_A+1
            elif x[i][j]=='T':
                count_T=count_T+1
            elif x[i][j]=='C':
                count_C=count_C+1
            elif x[i][j]=='G':
                count_G=count_G+1
            else: 
                count_other=count_other+1
        pos = j
        total=count_A+count_T+count_C+count_G+count_other
        freq_a=float(count_A)/float(total)
        freq_t=float(count_T)/float(total)
        freq_c=float(count_C)/float(total)
        freq_g=float(count_G)/float(total)
        freq_other=float(count_other)/float(total)
        result_a = (pos, 'a', freq_a)
        result_c = (pos, 'c', freq_c)
        result_g = (pos, 'g', freq_g)
        result_t = (pos, 't', freq_t)
        results=(result_a, result_c, result_g, result_t)
        [position.append(result) for result in results]
    return(position) 

def output_csv(outname,y):
    """ Function to output results to CSV """
    header = ['pos','seq','freq']
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

#x = fastq_reader('227_CCGTCC_L006_R1_002.fastq')
#y = position(x)
#output_csv('227_CCGTCC_L006_R1_002_pos_bias.csv',y)
