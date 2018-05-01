# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 11:25:39 2017

@author: Lucas Boatwright
"""

from sys import argv, stderr, exit
import argparse
from datetime import datetime
import subprocess as sp
from collections import OrderedDict
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio import SeqIO

def parse_arguments():
    """Parse arguments passed to script"""
    parser = argparse.ArgumentParser(description=
    "This script is designed to identify Common Orthologous REgions (COREs) \
    \ngiven either two parsed and filtered BLAST files and the FASTA files \
    \ncontaining the corresponding sequences (see 'rbh' subcommand) or given a \
    \nCSV with orthologous pairs and the FASTA files containing the \
    \ncorresponding sequences (see 'known' subcommand). If 'rbh' is used, \
    \nBLAST files are assumed to have a single best hit for each query. \
    \nCOREs are found in Species<AB>-species<BA>_ortholog_COREs.bed.\n\
    \n\tExample (known): python {0} -o orthologs.csv -f A.fasta -g B.fasta\
    \n\tExample (rbh): python {0} -a parsed_BLAST_A -b parsed_BLAST_B -f A.fasta -g B.fasta\
    \n\
    \nExample BLAST_A: \
    \n\tcontig_10340	101674	contig_252388	66445	1621	0.	80.4	80.4878048780488	533	32644	32124	65919	66445	-\
    \n\nExample BLAST_B: \
    \n\tcontig_252388	66865	contig_10340	101674	280762	0.	91.6	91.6699410609037	66170	66192	58	32362	98486	-\
    \n".format(argv[0]),
    formatter_class = argparse.RawDescriptionHelpFormatter)
    
    subparsers = parser.add_subparsers(help='help for subcommand')
    
    # Parser for known orthologous pairs, if available
    parse_known = subparsers.add_parser('known', help='If orthologous pairs \
    are known, use a user-defined CSV with Species_A_GeneID,Species_B_GeneID')

    required_known = parse_known.add_argument_group('required arguments')   
    
    required_known.add_argument('-o', "--ORTHOLOGS", type=str, required=True,
                             help='CSV file: Species_A_GeneID,Species_B_GeneID')
                             
    required_known.add_argument("-f", "--FASTA_A", type=str, required=True,
    help="FASTA file of all sequences for species A.", 
    action="store")
    
    required_known.add_argument("-g", "--FASTA_B", type=str, required=True,
    help="FASTA file of all sequences for species B.", 
    action="store")

    # Parser for performing reciprocal best-hit approach if orthologs unknown
    parse_rbh = subparsers.add_parser("rbh", help="If orthologs are not yet \
    identified, use a reciprocal best-hit approach to determine orthologous \
    pairs.")
    
    required_rbh = parse_rbh.add_argument_group('required arguments')    
    
    required_rbh.add_argument("-a", "--BLAST_A", type=str, required=True,
    help="Parsed BLAST file from species A in tab-separated format.", 
    action="store")
    
    required_rbh.add_argument("-b", "--BLAST_B", type=str, required=True,
    help="Parsed BLAST file from species B in tab-separated format.", 
    action="store")
    
    required_rbh.add_argument("-f", "--FASTA_A", type=str, required=True,
    help="FASTA file of all sequences for species A.", 
    action="store")
    
    required_rbh.add_argument("-g", "--FASTA_B", type=str, required=True,
    help="FASTA file of all sequences for species B.", 
    action="store")

    return parser.parse_args() 

#==== Multiple Sequence Alignment Functions ====#
def get_blast_dict(blast):
    """Generate a dictionary for a BLAST file."""
    blast_dict = {}
    with open(blast) as file:
        for line in file:
            split_line = line.split()
            blast_dict[split_line[0]] = split_line[2]
    return blast_dict


def reciprocated_contigs(blast_a_dict, blast_b_dict):
    """Determine which contigs hit reciprocally."""
    for key, value in blast_a_dict.items():
        try:
            if (key == blast_b_dict[value]):
                yield (key, value)
        except KeyError:
            stderr.write("\tNo key for {0} in BLAST_B.\n".format(value))
       
            
def perform_alignment():
    """Run MAFFT to obtain a multiple sequence alignment between orthologs"""
    sp.call(['/apps/mafft/7.127/bin/mafft', '--adjustdirection', '--clustalout', '--preservecase',
             'ortho_seqs.fasta'], 
             stdout = open("mafft_temp.out",'w'), 
             stderr = open("mafft_temp.log",'w'))
    sp.call(['cat','mafft_temp.out'],
            stdout = open("all_mafft_alignments.txt",'a'))
    sp.call(['cat','mafft_temp.log'], 
            stdout = open("all_mafft_logs.log",'a'))


def obtain_consensus(pair):
    """Obtain a consensus from the multiple sequence alignment"""
    consensus_sequences = open("consensus_sequences.fasta","a")
    align=AlignIO.read("mafft_temp.out","clustal")
    summary_align = AlignInfo.SummaryInfo(align)
    consensus = summary_align.dumb_consensus(threshold=0.51, ambiguous='N')
    new_sequence_name = ">{0},{1}\n{2}\n".format(pair[0], pair[1], consensus)
    consensus_sequences.write(new_sequence_name)  
    
    
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
    
    
def fasta_to_lists(fasta_file):
    """Isolate deflines and sequences from FASTA"""
    with open(fasta_file, "r") as f:
        file = f.readlines()
    deflines = []
    sequences = []
    sequence = ""
    for line in file:
        if line.startswith(">"):
            deflines.append(line.rstrip().lstrip('>'))
            if sequence:
                sequences.append(sequence)
                sequence = ""
        else:
            sequence += line.rstrip()
    sequences.append(sequence)
    return deflines, sequences
    
    
def read_ortholog_csv(orthologs):
    """If given a table of orthologs in CSV, yield pairs for alignment"""
    with open(orthologs) as file:
        for line in file:
            pair = [i.rstrip() for i in line.split(',')]
            yield pair
    
    
def write_and_align_orthologs(fasta_a_dict, fasta_b_dict, 
                              reciprocated_contigs):
    """Write ortholog pairs to file and run MAFFT on pair"""
    species_a_orthologs = open("SpeciesA-speciesB_orthologs.fasta",'w')
    species_b_orthologs = open("SpeciesB-speciesA_orthologs.fasta",'w')
    pairs = open("reciprocated_blast_hits.txt",'w')    
    for pair in reciprocated_contigs:
        if pair:
            pairs.write("{0},{1}\n".format(pair[0],pair[1]))
            with open("ortho_seqs.fasta",'w') as mafft_input:
                mafft_input.write(">{0}\n{1}\n>{2}\n{3}\n".format(pair[0],
                                                      fasta_a_dict[pair[0]],
                                                      pair[1],
                                                      fasta_b_dict[pair[1]]))
            species_a_orthologs.write(">{0}\n{1}\n".format(pair[0],
                                                      fasta_a_dict[pair[0]])) 
            species_b_orthologs.write(">{0}\n{1}\n".format(pair[1],
                                                      fasta_b_dict[pair[1]]))                                          
            perform_alignment()
            obtain_consensus(pair)
    pairs.close()
    species_a_orthologs.close()
    species_b_orthologs.close()
    
    
def alignment_to_fasta():
    """Convert multiple alignment file --> FAA to FASTA. Alignment was CLUSTAL 
    output format for visualization purposes."""
    open("all_mafft_alignments.faa",'w').close()
    alignments = list(AlignIO.parse("all_mafft_alignments.txt", "clustal"))

    for i in range(0,len(alignments)):
        with open("all_mafft_alignments.faa", "a") as handle:
            SeqIO.write(alignments[i], handle, "fasta")    
            
#==== BED GENERATION FUNCTIONS ====#
def get_start_index(msa_seqs):
    """Get start indices from MSA"""
    start_indices = []
    for sequence in msa_seqs:
        for nucleotide in sequence:
            if nucleotide != '-':
                start_indices.append(sequence.index(nucleotide))
                break
    return start_indices


def get_stop_index(msa_seqs):
    """Get stop indices from MSA"""
    stop_indices = []
    for sequence in msa_seqs:
        reverse_seq = sequence[::-1]
        for nucleotide in reverse_seq:
            if nucleotide != '-':
                index=reverse_seq.index(nucleotide)
                stop = len(reverse_seq) - index
                stop_indices.append(stop)
                break
    return stop_indices
    
    
def get_strand(msa_deflines):
    """Get contig strand (+/-) from MSA deflines"""
    strands = []
    for definition in msa_deflines:
        if "_R_" not in definition:
            strands.append("+")
        elif "_R_" in definition:
            strands.append("-")
    return strands
    
    
def determine_segment_distribution(ortholog_one, ortholog_two):
    """Determine which contig occurs first"""
    segment_number = 1       
    if ortholog_one[1] > ortholog_two[1]:
        #species two starts first
        line = "{0}\t{1}\t{2}\t{0}|{4}|{3}".format(ortholog_one[3], \
        ortholog_two[1], ortholog_one[1], ortholog_two[0], segment_number)
        segment_number += 1
        if ortholog_one[2] > ortholog_two[2]:
            #species two stops first
            line = "{0}\n{3}\t{1}\t{2}\t{3}|{4}|CORE".format(line, \
            ortholog_one[1], ortholog_two[2], ortholog_one[3], segment_number)
            segment_number += 1
            line = "{0}\n{3}\t{1}\t{2}\t{3}|{5}|{4}".format(line, \
            ortholog_two[2],ortholog_one[2], ortholog_one[3], ortholog_one[0], segment_number)
        elif ortholog_one[2] < ortholog_two[2]:
            #species one stops first
            line = "{0}\n{3}\t{1}\t{2}\t{3}|{4}|CORE".format(line, \
            ortholog_one[1], ortholog_one[2], ortholog_one[3], segment_number)
            segment_number += 1
            line = "{0}\n{3}\t{1}\t{2}\t{3}|{5}|{4}".format(line, \
            ortholog_one[2],ortholog_two[2], ortholog_one[3], ortholog_two[0], segment_number)
        elif ortholog_one[2] == ortholog_two[2]:
            #species are the exact same length
            line = "{0}\n{3}\t{1}\t{2}\t{3}|{4}|CORE".format(line, \
            ortholog_one[1], ortholog_one[2], ortholog_one[3], segment_number)
            #line = "{0}\n{2}\t{1}\t{1}\t{2}_empty".format(line, \
            #ortholog_one[2], ortholog_one[3])
    elif ortholog_one[1] < ortholog_two[1]:
        #species one starts first
        line = "{0}\t{1}\t{2}\t{0}|{4}|{3}".format(ortholog_one[3], \
        ortholog_one[1], ortholog_two[1], ortholog_one[0], segment_number)
        segment_number += 1
        if ortholog_one[2] > ortholog_two[2]:
            #species two stops first
            line = "{0}\n{3}\t{1}\t{2}\t{3}|{4}|CORE".format(line, \
            ortholog_two[1], ortholog_two[2], ortholog_one[3], segment_number)
            segment_number += 1
            line = "{0}\n{3}\t{1}\t{2}\t{3}|{5}|{4}".format(line, \
            ortholog_two[2],ortholog_one[2], ortholog_one[3], ortholog_one[0], segment_number)
        elif ortholog_one[2] < ortholog_two[2]:
            #species one stops first
            line = "{0}\n{3}\t{1}\t{2}\t{3}|{4}|CORE".format(line, \
            ortholog_two[1], ortholog_one[2], ortholog_one[3], segment_number)
            segment_number += 1
            line = "{0}\n{3}\t{1}\t{2}\t{3}|{5}|{4}".format(line, \
            ortholog_one[2],ortholog_two[2], ortholog_one[3], ortholog_two[0], segment_number)
        elif ortholog_one[2] == ortholog_two[2]:
            #species are the exact same length
            line = "{0}\n{3}\t{1}\t{2}\t{3}|{4}|CORE".format(line, \
            ortholog_two[1], ortholog_two[2], ortholog_one[3], segment_number)
            #line = "{0}\n{2}\t{1}\t{1}\t{2}_empty".format(line, \
            #ortholog_one[2], ortholog_one[3])       
    elif ortholog_one[1]==ortholog_two[1]:
        #same start site and line should be empty
        #line = "{0}\t{1}\t{2}\t{0}_empty".format(ortholog_one[3], \
        #ortholog_one[1], ortholog_two[1])
        segment_number += 1
        if ortholog_one[2] > ortholog_two[2]:
            #species two stops first
            line = "{3}\t{1}\t{2}\t{3}|{0}|CORE".format(segment_number, \
            ortholog_two[1], ortholog_two[2], ortholog_one[3])
            segment_number += 1
            line = line + "\n{3}\t{1}\t{2}\t{3}|{0}|{4}\n".format(segment_number, \
            ortholog_two[2],ortholog_one[2], ortholog_one[3], ortholog_one[0])
        elif ortholog_one[2] < ortholog_two[2]:
            #species one stops first
            line = "{3}\t{1}\t{2}\t{3}|{0}|CORE".format(segment_number, \
            ortholog_two[1], ortholog_one[2], ortholog_one[3])
            segment_number += 1
            line = line + "\n{3}\t{1}\t{2}\t{3}|{0}|{4}\n".format(segment_number, \
            ortholog_one[2],ortholog_two[2], ortholog_one[3], ortholog_two[0])
        elif ortholog_one[2] == ortholog_two[2]:
            #species are the exact same length
            line = "{3}\t{1}\t{2}\t{3}|{0}|CORE".format(segment_number, \
            ortholog_one[1], ortholog_one[2], ortholog_one[3])
            #line = "{0}\n{2}\t{1}\t{1}\t{2}_empty".format(line, \
            #ortholog_one[2], ortholog_one[3])        
    return line


def pair_ortholog_information(fasta_deflines, msa_deflines, start_indices, \
stop_indices):
    """Group contig name with start and stop"""
    stop = len(msa_deflines) - 1
    first = 0
    second = 1
    fasta_deflines_index=0
    
    line_list = []
    
    while second <= stop:
        # take a pair at a time for BED file line generation
        if "_R_" not in msa_deflines[first]:
            ortholog_one = [msa_deflines[first].split("_")[0].strip(">"), \
            start_indices[first], stop_indices[first], \
            fasta_deflines[fasta_deflines_index].strip(">")]
        else:
            species = msa_deflines[first].lstrip(">_R_").split("_")[0]
            ortholog_one = [species, \
            start_indices[first], stop_indices[first], \
            fasta_deflines[fasta_deflines_index].strip(">")]            
        
        if "_R_" not in msa_deflines[second]:
            ortholog_two = [msa_deflines[second].split("_")[0].strip(">"), \
            start_indices[second], stop_indices[second]]
        else:
            species = msa_deflines[second].lstrip(">_R_").split("_")[0]
            ortholog_two = [species, \
            start_indices[second], stop_indices[second]]
        
        line = determine_segment_distribution(ortholog_one, ortholog_two)
        line_list.append(line)
        
        first += 2
        second += 2
        fasta_deflines_index += 1

    return line_list


def write_bed_file(line_list):
    """Write BED file"""
    output = open("ortholog_divisions.bed", "w")
    for line in line_list:
        output.write(line.rstrip() + "\n")

            
def faa_to_bed():
    """Generate BED files from the multiple sequence alignment output.""" 
    fasta_deflines, fasta_seqs = fasta_to_lists("consensus_sequences.fasta")
    msa_deflines, msa_seqs = fasta_to_lists("all_mafft_alignments.faa")
    start_indices = get_start_index(msa_seqs)
    stop_indices = get_stop_index(msa_seqs)
    line_list = pair_ortholog_information(fasta_deflines, msa_deflines, 
                                          start_indices, stop_indices)                                          
    write_bed_file(line_list)

    sp.call(["sed",'s/,/|/g', "ortholog_divisions.bed"], 
            stdout = open("temp.bed",'w'))
    sp.call(["mv","temp.bed","ortholog_divisions.bed"])
       
       
# Functions for species-specific BED files from COREs    
def open_bed_file(bed_file_name):
    """Open BED file of sequence overlaps"""
    with open(bed_file_name) as f:
        cores = [i.split() for i in f.readlines() if "CORE" in i]
    return cores
    
    
def sort_ortholog_information(faa_defs, faa_seqs):
    """Group contig name with start and stop"""
    first_species_defs = []
    second_species_defs = []
    first_species_seqs = []
    second_species_seqs = []
    for x, i in enumerate(faa_defs):
        if x%2==0:
            first_species_defs.append(i)
            first_species_seqs.append(faa_seqs[x])
        elif x%2==1:
            second_species_defs.append(i)
            second_species_seqs.append(faa_seqs[x])
    # species_defs can only be used for their strand information
    return first_species_defs, first_species_seqs, second_species_defs, second_species_seqs
    
    
def reverse_complement(seq):
    """Get reverse complement of a given sequence"""
    comp_dict = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
    return "".join([comp_dict[base] for base in reversed(seq)]) 


def get_overlaps(cores, first_seqs, sec_seqs):
    """Used BED overlap coords to isolate the corresponding seqs"""
    first_dict = OrderedDict()
    second_dict = OrderedDict()
    for x,i in enumerate(cores):
        first_dict[i[0].split("|")[0]]  = first_seqs[x][int(i[1]):int(i[2])]
        second_dict[i[0].split("|")[1]] = sec_seqs[x][int(i[1]):int(i[2])]
    return first_dict, second_dict
    
    
def remove_gaps(species_overlaps):
    """Remove '-' from all sequences"""
    corrected_overlap = OrderedDict()
    for key, value in species_overlaps.items():
        try:
            corrected_overlap[key]= "".join(value.replace('-',''))
        except ValueError:
            corrected_overlap[key]=value
    # Called 'corrected_overlap' since it may have undergone 
    # reverse complementation via MAFFT
    return corrected_overlap
    

def get_coords(gapless_overlaps, species_fasta, species_strands):
    """Using the isolated overlaps (dict), perform reverse comp where necessary 
    and obtain the coords of the original overlap in its full-length sequence"""
    fasta = OrderedDict()
    # species_fasta represents the original, full-length sequences
    defs, seqs = fasta_to_lists(species_fasta)
    for x, i in enumerate(defs):
        fasta[i]=seqs[x]
    totals = []
    indices = []
    count = 0
    # need to figure out species dict, fasta, strand information
    for key, value in gapless_overlaps.items():
        strand = species_strands[count]
        if strand=='-':
            indices.append( (strand, 
            seqs[count].index(reverse_complement(value) ), 
            seqs[count].index(reverse_complement(value)) + len(value)) )
            totals.append( fasta[key].count(reverse_complement(value)) )
        elif strand=='+':
            indices.append( (strand, seqs[count].index(value), 
            seqs[count].index(value) + len(value) ) )
            totals.append( fasta[key].count(value) )
        count+=1
    try: assert totals.count(1)==len(gapless_overlaps)
    except AssertionError:
        print(totals.count(1))
        print(len(gapless_overlaps))
        exit(1)
    assert len(gapless_overlaps)==len(indices)
    assert len(indices)==len(species_strands)
    # indices has strand, start and stop information
    return indices, defs
    

def write_species_beds(indices, fasta_name, overlaps, species_defs):
    """Write bed file with coords, strands, seq names"""
    # indices = list of tuples
    # fasta_name = string
    # overlaps = list of lists
    # species_defs = list of strings
    
    # WRT = 'with respect to'
    output_name = fasta_name.split('_')[0] + "_ortholog_COREs.bed"
    output = open(output_name,'w')
    assert len(indices)==len(overlaps)==len(species_defs)
    for x, i in enumerate(indices):
        output.write("{0}\t{1}\t{2}\t{3}\t0\t{4}\n".format(
        species_defs[x], i[1], i[2], overlaps[x][0], i[0]))
        
    
def get_species_bed_files():
    """Generate a BED file for each species from the COREs"""
    # For any confusion that may ensue from brevity,
    # seqs is shorthand for sequences
    # defs = FASTA definition lines
    # coords = coordinates
    fasta_a = "SpeciesA-speciesB_orthologs.fasta"
    fasta_b = "SpeciesB-speciesA_orthologs.fasta"
    cores = open_bed_file("ortholog_divisions.bed")
    faa_defs, faa_seqs = fasta_to_lists("all_mafft_alignments.faa")
    first_species_defs, first_species_seqs, second_species_defs, \
        second_species_seqs = sort_ortholog_information(
                                                    faa_defs, faa_seqs)
    first_overlaps, second_overlaps = get_overlaps(
                        cores, first_species_seqs, second_species_seqs)
                        
    first_gapless_overlaps  = remove_gaps(first_overlaps)
    second_gapless_overlaps = remove_gaps(second_overlaps)    

    first_strands = get_strand(first_species_defs)
    second_strands= get_strand(second_species_defs)        

    first_indices, first_original_defs = get_coords(
        first_gapless_overlaps, fasta_a, first_strands)
    second_indices, second_original_defs = get_coords(
    second_gapless_overlaps, fasta_b, second_strands)
    # write the species BED files
    write_species_beds(first_indices, fasta_a, cores, first_original_defs)
    write_species_beds(second_indices, fasta_b, cores, second_original_defs)                                        
    
    
if __name__ == "__main__":
    start = datetime.now()
    args = parse_arguments()
    try:
        stderr.write("Executed: python {0} -a {1} -b {2} -f {3} -g {4}\n".format(
            argv[0],
            args.BLAST_A,
            args.BLAST_B,
            args.FASTA_A,
            args.FASTA_B))
    
        stderr.write("Generating BLAST_A dictionary.\n")
        blast_a_dict = get_blast_dict(args.BLAST_A)    
        
        stderr.write("Generating BLAST_B dictionary.\n")
        blast_b_dict = get_blast_dict(args.BLAST_B)

        stderr.write("Generating reciprocal best-hit pairs.\n")
        recip_contigs = reciprocated_contigs(blast_a_dict, blast_b_dict)
    
    except AttributeError:
        stderr.write("Executed: python {0} -o {1} -f {2} -g {3}\n".format(
            argv[0],
            args.ORTHOLOGS,
            args.FASTA_A,
            args.FASTA_B))
            
        stderr.write("Reading orthologs from CSV.\n")
        recip_contigs = read_ortholog_csv(args.ORTHOLOGS)
    
    stderr.write("Generating FASTA_A dictionary.\n")
    fasta_a_dict = fasta_to_dict(args.FASTA_A)
    
    stderr.write("Generating FASTA_B dictionary.\n")
    fasta_b_dict = fasta_to_dict(args.FASTA_B)
    
    stderr.write("Performing sequence alignment and writing tables.\n" + 
                "\t-- This may take a while.\n")
    write_and_align_orthologs(fasta_a_dict, fasta_b_dict,
                              recip_contigs)
                   
    stderr.write("Converting CLUSTAL output to FAA.\n")               
    alignment_to_fasta()
    
    stderr.write("Generating CORE BED file from MSA.\n")
    faa_to_bed()
    
    stderr.write("Generating Species-specific BED files from CORE BED.\n")
    get_species_bed_files()
    
    stop = datetime.now()
    stderr.write("Runtime: {0}\n".format(str(stop - start)))
