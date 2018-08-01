# -*- coding: utf-8 -*-
"""
Created on Fri May  5 11:30:16 2017

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
    "This script is designed to take two bed files and generate a mean-difference plot of the CORE lengths\n\n \
    Example: python {0} -a A.bed -b B.bed -o BED_mean_difference.pdf".format(argv[0]),
    formatter_class = argparse.RawDescriptionHelpFormatter)
    
    requiredNamed = parser.add_argument_group('required arguments')

    requiredNamed.add_argument("-a", "--BED_A", type=str, required=True,
    help="BED file A", action="store")

    requiredNamed.add_argument("-b", "--BED_B", type=str, required=True,
    help="BED file B", action="store")

    parser.add_argument("-o", "--OUTPUT", type=str, required=False,
    help="Output file name.", action="store")

    return parser.parse_args() 
        

def mean_diff_cores(bed_a_filename, bed_b_filename, 
                    output="BED_mean_difference.pdf"):
    bed_a = pd.read_table(bed_a_filename, names=["seq","start","stop","commonID","fill","strand"])
    bed_b = pd.read_table(bed_b_filename, names=["seq","start","stop","commonID","fill","strand"])

    bed_a["bed_a_cores"] = bed_a["stop"] - bed_a["start"]
    bed_b["bed_b_cores"] = bed_b["stop"] - bed_b["start"]

    bed_a_cores = pd.DataFrame(bed_a["bed_a_cores"])
    bed_a_cores["bed_b_cores"] = bed_b["bed_b_cores"]

    bed_a_cores["Mean"]=bed_a_cores.mean(axis=1)
    bed_a_cores["Difference"]=bed_a_cores["bed_a_cores"] - bed_a_cores["bed_b_cores"]
    bed_a_cores["commonID"] = bed_a["commonID"]

    print(bed_a_cores["Difference"].describe())
    print(bed_a_cores["Mean"].describe())

    high_diff = bed_a_cores[bed_a_cores["Difference"] >= 1000]
    print(high_diff["commonID"])

    plt.figure()
    plt.scatter(bed_a_cores["Mean"],bed_a_cores["Difference"], facecolors='none')
    plt.xlabel("Length Mean (bp)")
    plt.ylabel("Length Difference (bp)")
    plt.savefig(output, format="pdf")
    plt.show()        
        
    
if __name__ == "__main__":
    start = datetime.now()
    args = parse_arguments()
    if args.OUTPUT:
        stderr.write("Executed: python {0} -a {1} -b {2} -o {3}\n".format(
                    argv[0],args.BED_A, args.BED_B, args.OUTPUT))
        mean_diff_cores(args.BED_A, args.BED_B, 
                        args.OUTPUT)
    else:
        stderr.write("Executed: python {0} -a {1} -b {2}\n".format(argv[0],
                 args.BED_A, args.BED_B))
        mean_diff_cores(args.BED_A, args.BED_B)
    stop = datetime.now()
    stderr.write("Runtime: {0}\n".format(str(stop - start)))

