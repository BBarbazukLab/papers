
__author__="Lucas Boatwright"
__date__ ="$Oct 3, 2014 3:07:06 PM$"

import sys

def fasta_to_lists(file_name):
    """Isolate deflines and sequences from FASTA"""
    with open(file_name, "r") as f:
        file = f.readlines()
    deflines = []
    sequences = []
    sequence = ""
    for line in file:
        if line.startswith(">"):
            deflines.append(line.rstrip().lstrip(">"))
            if sequence:
                sequences.append(sequence)
                sequence = ""
        else:
            sequence += line.rstrip()
    sequences.append(sequence)
    return deflines, sequences

def generate_table_of_changes_and_new_fasta(file_name, deflines, seqs, new_prefix):
    """Write table of old and new sequence names"""
    old_name = file_name.split(".fa")[0]
    new_fasta_name = old_name + ".renamed.fasta"
    table_name = old_name + ".renamed.table"
    new_fasta = open(new_fasta_name, "w")
    table = open(table_name, "w")
    table.write(file_name + "\t" + new_fasta_name + "\n")
    for x,i in enumerate(deflines):
        table.write(i + "\t" + new_prefix + "_" + str(x+1) + "\n")
        new_fasta.write(">" + new_prefix + "_" + str(x+1) + "\n")
        new_fasta.write(seqs[x] + "\n")

def main():
    """Execute all functions"""
    try:
        deflines, seqs = fasta_to_lists(sys.argv[1])
        generate_table_of_changes_and_new_fasta(sys.argv[1], deflines, seqs, sys.argv[2])
    except IndexError:
        print "\n\npython rename_fasta_sequences_and_log_changes.py old_fasta_file.fa new_prefix\n\n"
        sys.exit(1)

if __name__ == "__main__":
    main()

