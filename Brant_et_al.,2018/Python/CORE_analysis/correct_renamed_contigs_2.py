import sys

if len(sys.argv) != 3:
    print("\npython correct_renamed_contigs_2.py gene_trans.txt Fasta.renamed.table\n")
    sys.exit(1)

with open(sys.argv[1]) as f:
    orthologs = f.readlines()

with open(sys.argv[2]) as f:
    table = f.readlines()

table_dict = {}
for line in table:
    split_line = line.split()
    table_dict[split_line[1].rstrip()]=split_line[0].rstrip()

output = open(sys.argv[1].split(".")[0] + "_corrected.txt","w")

for line in orthologs:
    split_line = line.split()
    output.write(split_line[0] + "\t" + table_dict[split_line[1]] + "\n")

output.close()

