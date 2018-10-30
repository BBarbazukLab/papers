"""genes_from_edge_table.py: This script pulls a single gene from the Coexpression edge table to identify adjacent nodes and prints the gene with neighbors to the screen."""

__author__      = "J. Lucas Boatwright"

import pandas as pd
from sys import argv, exit

if len(argv) != 3:
    print("\npython {0} edge.table gene_symbol\n".format(argv[0]))
    exit(1)

df = pd.read_csv(argv[1])
df["name"] = df["name"].apply(lambda x: x.lower())
df = df[df["name"].str.contains(argv[2].lower())]
df["first"] = df["name"].apply(lambda x: x.split(" (interacts with) ")[0])
df["second"] = df["name"].apply(lambda x: x.split(" (interacts with) ")[1])
df = df.loc[(df["first"]==argv[2].lower()) | (df["second"]==argv[2].lower())]

genes = pd.concat([df["first"], df["second"]])
genes = set(genes)
for gene in genes:
    print(gene)
