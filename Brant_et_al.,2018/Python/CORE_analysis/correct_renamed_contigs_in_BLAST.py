import sys

def usage():
    if len(sys.argv) != 4:
        print("\npython correct_renamed_contigs_in_BLAST.py Blast_output name_conversion.table Trinity.renamed.table\n")
        sys.exit(1)

def open_blast(file_name):
    with open(file_name) as f:
        orthologs = f.readlines()
    return orthologs


def get_table_dict(file_name):
    with open(file_name) as f:
        table = f.readlines()

    table_dict = {}
    for line in table:
        split_line = line.split('"')
        table_dict[split_line[0].rstrip()]=split_line[1].rstrip()
    return table_dict


def get_trinity_dict(file_name):
    """Generate a dict from Trinity.renamed.table"""
    with open(file_name) as f:
        table = f.readlines()

    table_dict = {}
    for line in table:
        split_line = line.split()
        table_dict[split_line[1].rstrip()]=split_line[0].rstrip()
    return table_dict


def replace_trinity_names(blast, trinity_table):
    new_blast = []
    for line in blast:
        split_line = line.split()
        new_blast.append("\t".join(split_line[0:2] + [trinity_table[split_line[2]]] + split_line[3:]))
    return new_blast


def write_output(blast_file_name, table_dict, orthologs):
    """Write output"""
    output = open(blast_file_name.split(".")[0] + "_corrected.txt","w")

    for line in orthologs:
#        split_line = line.split()
#        try:
#            out = table_dict[split_line[0]] + "\t" + "\t".join(split_line[1:]) + "\n"
#        except KeyError:
#            continue
        output.write(line + "\n")
    output.close()


if __name__ == "__main__":
    usage()
    blast = open_blast(sys.argv[1])
    gene_sym_table = get_table_dict(sys.argv[2])
    trinity_table = get_trinity_dict(sys.argv[3])
    new_blast = replace_trinity_names(blast, trinity_table)
    write_output(sys.argv[1], gene_sym_table, new_blast)
