# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 13:37:06 2017

@author: Lucas Boatwright
"""

from sys import argv, stderr
import argparse
from datetime import datetime
import subprocess as sp
from itertools import chain
from collections import OrderedDict

def parse_arguments():
    """Parse arguments passed to script"""
    parser = argparse.ArgumentParser(description=
    "This script takes each group from grouped_hits.txt and executes \
    \nCAP3 on each group in an attempt to merge transcripts and remove \
    \nredundant contigs. Requires a FASTA file containing all sequences from \
    \ngrouped_hits.txt. Auto-detects 'grouped_hits.txt' in the current \
    \ndirectory. All contigs after compression are in 'contigs_collapsed.fasta'\n\n \
    Example: python {0} -d contigs.fa".format(argv[0]),
    formatter_class = argparse.RawDescriptionHelpFormatter)

    requiredNamed = parser.add_argument_group('required arguments')

    requiredNamed.add_argument("-f", "--FASTA", type=str, required=True,
    help="FASTA file containing all sequences within grouped_hits.txt.", 
    action="store")

    return parser.parse_args() 
    
    
def fasta_to_dict(fasta_file):
    """Consolidate deflines and sequences from FASTA as dictionary"""
    deflines = []
    sequences = []
    sequence = ""
    with open(fasta_file, "r") as file:
        for line in file:
                if line.startswith(">"):
                        deflines.append(line.rstrip().lstrip('>'))
                        if sequence:
                                sequences.append(sequence)
                                sequence = ""
                else:
                        sequence += line.rstrip()
        sequences.append(sequence)
    fasta_dict = OrderedDict()
    for x, defline in enumerate(deflines):
        fasta_dict[defline]=sequences[x]
    return fasta_dict
        

def contig_compression_looper(fasta_dict):
    """Attempt to collapse groups of isoforms with CAP3"""	
    with open("grouped_hits.txt","r") as f:	
        file = f.readlines()
		
    # output file for initial contig #s and output numbers (this will require renaming)	
    contig_compression = open("contig_compression_data.txt","w")
    contig_in_out = open("contigs_in_out.txt","w")
    new_contigs = open("contigs_collapsed.fasta","w")

    # for renaming contigs
    group_num = 0
    for i in file:
        contigs_in = i.rsplit()
        pre_counts = 0
        group_num+=1
        group_deci = 1
        open("cap3_temp",'w').close()
  
        for j in contigs_in:
            group = j.rsplit("_")			
            group.pop(-1)
            contig_in_out.write("{0}\t".format(j))
            contig_in_out.flush()
            with open("cap3_temp",'a') as output:
                output.write(">{0}\n{1}\n".format(j, fasta_dict[j]))
            pre_counts+=1
            
        sp.call(['/apps/cap3/20120705/bin/cap3', 'cap3_temp', 
                 '-o', '25', '-p', '80'], 
                stdout=open("cap3_output",'w'), 
                stderr=open("cap3_log",'w'))
	
        with open("cap3_temp.cap.contigs","r") as f:
            contig_counts = f.readlines()
			
        c_counts = 0
        s_counts = 0	
        new_defline = "_".join(group)
        contig_in_out.write("||\t")
			
        for i in contig_counts:
            if ">" in i:
                c_counts+=1
                new_contigs.write(">{0}_{1}.{2}\n".format(new_defline,
                                                          group_num,
                                                          group_deci))
                contig_in_out.write("{0}_{1}.{2}\t".format(new_defline,
                                                            group_num,
                                                            group_deci))
                group_deci += 1
            else:
                new_contigs.write(i)
		
        contig_compression.write("pre_cap3: {0}\tcontigs: {1}".format(pre_counts,c_counts))
        contig_compression.flush()
		
        with open("cap3_temp.cap.singlets","r") as f:
            singlet_counts = f.readlines()
			
        for i in singlet_counts:
            if ">" in i:
                s_counts+=1
                new_contigs.write(">{0}_{1}.{2}\n".format(new_defline,group_num,group_deci))
                contig_in_out.write("{0}_{1}.{2}\t".format(new_defline,group_num,group_deci))				
                group_deci += 1
            else:
                new_contigs.write(i)				
	
        total_post = c_counts + s_counts
	
        contig_compression.write("\tsinglets: {0}\ttotal_post_cap3: {1}\tdifference: {2}"\
        .format(s_counts,total_post,(pre_counts-total_post)))
		
        contig_compression.flush()
        contig_in_out.write("\n")
        contig_in_out.flush()	
        new_contigs.flush()
        contig_compression.write("\n")
        
        sp.call(['cat','cap3_temp.cap.ace'], stdout=open("ace_files.ace",'a'))
        sp.call(['cat','cap3_output'], stdout=open("compressed_cap3_output.txt",'a'))
        sp.call(['cat','cap3_log'], stdout=open("compressed_cap3_log.txt",'a'))
        
    sp.call(['rm','cap3_temp.cap.singlets','cap3_temp.cap.info',
             'cap3_temp.cap.contigs','cap3_temp.cap.contigs.qual',
             'cap3_temp.cap.contigs.links','cap3_temp.cap.ace',
             'cap3_log','cap3_temp','cap3_output'])


def isolate_ungrouped_contigs(fasta_dict):
    """Determine which transcripts were not grouped and concatenate both
    contigs_collapsed.fasta and ungrouped contigs into new transcriptome"""
    with open("grouped_hits.txt") as f:
        grouped_contigs = list(chain.from_iterable(
                                [i.split() for i in f.readlines()]
                                                    )
                                )
    db_deflines = fasta_dict.keys()
    ungrouped_contigs = [i for i in db_deflines if i not in grouped_contigs]
    with open("ungrouped_contigs.fasta",'w') as output:
        for contig in ungrouped_contigs:
            output.write(">{0}\n{1}\n".format(contig, fasta_dict[contig]))   
    sp.call(['cat','contigs_collapsed.fasta','ungrouped_contigs.fasta'], 
            stdout=open("full_assembly_after_collapse.fasta",'w'))
    
    
if __name__ == "__main__":
    start = datetime.now()
    args = parse_arguments()
    stderr.write("Executed: python {0} -d {1}\n".format(argv[0],args.FASTA))
    
    stderr.write("Reading FASTA file to dictionary.\n")
    fasta_dict = fasta_to_dict(args.FASTA)
    
    stderr.write("Starting compression loops.\n")
    contig_compression_looper(fasta_dict)

    stderr.write("Finding ungrouped contigs and finalizing compressed " + 
                    "assembly.\n")    
    isolate_ungrouped_contigs(fasta_dict)

    stop = datetime.now()
    stderr.write("Runtime: {0}\n".format(str(stop - start)))
