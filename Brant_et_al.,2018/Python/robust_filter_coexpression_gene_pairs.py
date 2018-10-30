# coding: utf-8

"""robust_filter_coexpression_gene_pairs.py: This script pulls a list of genes from the Coexpression edge table where neighbors must also be in the list and prints the gene with corresponding neighbors to the screen."""

__author__      = "J. Lucas Boatwright"

import pandas as pd
import sys

if len(sys.argv) != 2:
    print("python {0} gene_list.txt".format(sys.argv[0]))
    print("Coexpression_edgeTable.csv MUST be in current directory.")
    sys.exit(1)

def pair_interactions(x):
    """Pair interacting nodes"""
    split_line = x.split()
    pair = "{0}_{1}".format(split_line[0], split_line[3])
    return pair


def read_genes(file):
    """Read in gene list from file"""
    with open(file) as f:
        genes = [i.rstrip().lower() for i in f.readlines()]
    return genes


if __name__ == "__main__":
    genes = read_genes(sys.argv[1])
    df = pd.read_csv("Coexpression_edgeTable.csv")
    df["shared name"] = df["shared name"].apply(lambda x: pair_interactions(x))
    df["A"] = df["shared name"].apply(lambda x: x.split("_")[0].lower())
    df["B"] = df["shared name"].apply(lambda x: x.split("_")[1].lower())

    filtered = df[df["A"].isin(genes) & df["B"].isin(genes)]
    output = sys.argv[1].split(".txt")[0]
    filtered.to_csv(output + ".csv", index=False)
