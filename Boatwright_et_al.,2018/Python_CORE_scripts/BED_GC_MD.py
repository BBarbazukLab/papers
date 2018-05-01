# -*- coding: utf-8 -*-
"""
Created on Tue May  9 09:39:17 2017

@author: Lucas Boatwright
"""

from sys import argv, stderr
import argparse
from datetime import datetime
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

def parse_arguments():
    """Parse arguments passed to script"""
    parser = argparse.ArgumentParser(description=
    "This script is designed to generate an MD plot of GC bias\n\n \
    Example: python {0} -a A.bed.stats.csv -b B.bed.stats.csv -o orthlogs.txt".format(argv[0]),
    formatter_class = argparse.RawDescriptionHelpFormatter)
    
    requiredNamed = parser.add_argument_group('required arguments')

    requiredNamed.add_argument("-a", "--BED_A_STATS", type=str, required=True,
    help="BED file A stats from pyfasta info", action="store")

    requiredNamed.add_argument("-b", "--BED_B_STATS", type=str, required=True,
    help="BED file B stats from pyfasta info", action="store")

    requiredNamed.add_argument("-o", "--ORTHOLOGS", type=str, required=True,
    help="Orthologs.", action="store")

    return parser.parse_args() 
        

def GC_MD_cores(bed_a_filename, bed_b_filename, orthologs):
    orthologs = pd.read_csv(orthologs, names = ["A","B"])
    
    bed_a = pd.read_csv(bed_a_filename, names=["A","A_GC"])
    bed_a["A"]=bed_a["A"].apply(lambda x: x.lstrip('>'))
    
    bed_b = pd.read_csv(bed_b_filename, names=["B","B_GC"])
    bed_b["B"]=bed_b["B"].apply(lambda x: x.lstrip('>'))
    
    all_df = pd.merge(orthologs, bed_a)
    all_df = pd.merge(all_df, bed_b)

    all_df["Mean"]=all_df[["A_GC","B_GC"]].mean(axis=1)
    all_df["Difference"]=all_df["A_GC"] - all_df["B_GC"]

    print(all_df["Difference"].describe())
    print(all_df["Mean"].describe())

    output="BED_GC_MD.pdf"
    plt.figure()
    plt.scatter(all_df["Mean"],all_df["Difference"], facecolors='none')
    plt.xlabel("%GC Mean")
    plt.ylabel("%GC Difference")
    plt.savefig(output, format="pdf")
    plt.show()        
        
    
if __name__ == "__main__":
    start = datetime.now()
    args = parse_arguments()
    stderr.write("Executed: python {0} -a {1} -b {2} -o {3}\n".format(
                    argv[0],args.BED_A_STATS, args.BED_B_STATS, args.ORTHOLOGS))
    GC_MD_cores(args.BED_A_STATS, args.BED_B_STATS, 
                        args.ORTHOLOGS)
    stop = datetime.now()
    stderr.write("Runtime: {0}\n".format(str(stop - start)))


