import sys

# Open BED file
with open(sys.argv[1]) as f:
    bed = f.readlines()

# Open Trinity.renamed.table
table_dict = {}
with open(sys.argv[2]) as f:
    for line in f:
        split_line = line.split()
        table_dict[split_line[1].rstrip()]=split_line[0].rstrip()

output = open(sys.argv[1].split(".bed")[0] + ".corrected.bed","w")

for line in bed:
    split_line = line.split()
    output.write(table_dict[split_line[0]] + "\t" + "\t".join(split_line[1:]) + "\n")

output.close()

