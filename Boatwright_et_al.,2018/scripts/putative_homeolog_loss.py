"""
Created on Tue May  9 13:48:20 2017

@author: Lucas Boatwright
"""

from sys import argv, stderr
import argparse
from datetime import datetime
import pandas as pd

def parse_arguments():
    """Parse arguments passed to script"""
    parser = argparse.ArgumentParser(description=
    "This script is designed to check for putative loss between homeologs \
    \nusing counts from homeolog COREs used in the PG analysis.\n\n \
    Example: python {0} -c ase_bayes_Cast2_lam_croc_flag.csv -b bayes_flag_sig_Cast2.csv".format(argv[0]),
    formatter_class = argparse.RawDescriptionHelpFormatter)
    
    requiredNamed = parser.add_argument_group('required arguments')

    requiredNamed.add_argument("-c", "--COUNTS", type=str, required=True,
    help="Bayes flagged CSV of read counts for Line and Tester", action="store")
    
    requiredNamed.add_argument("-b", "--BIAS", type=str, required=True,
    help="Bayes flagged CSV of loci demonstrating HSE", action="store")

    return parser.parse_args() 
    

def find_putative_loss(file_name):
    """Given file name, determine loci that are putatively lost based upon lack of expression"""
    df = pd.read_csv(file_name)
    output_base = file_name.split(".csv")[0]
    df["Line_mean"] = df[["LINE_TOTAL_1","LINE_TOTAL_2","LINE_TOTAL_3"]
                            ].mean(axis=1)
    df["Line_std"] = df[["LINE_TOTAL_1","LINE_TOTAL_2","LINE_TOTAL_3"]
                            ].std(axis=1)
    df["Tester_mean"] = df[["TESTER_TOTAL_1","TESTER_TOTAL_2","TESTER_TOTAL_3"]
                            ].mean(axis=1)
    df["Tester_std"]=df[["TESTER_TOTAL_1","TESTER_TOTAL_2","TESTER_TOTAL_3"]
                            ].std(axis=1)
    
    df["Line_zero"] = ((df["Tester_mean"] > 0 ) & (df["Line_mean"] == 0 ))
    df["Tester_zero"] = ((df["Line_mean"] > 0 ) & (df["Tester_mean"] == 0 ))      

    return df

    
def verify_HSE(df, hse_file):
    """Verify that each locus does not overlap 0.5 based upon 95% HDI"""
    flag_sig = pd.read_csv(hse_file)
    sig_true = flag_sig.columns[-1]
    merged_df = pd.merge(df, flag_sig)    
    output_df = merged_df[["commonID","Line_zero","Tester_zero","Line_mean","Tester_mean",sig_true]]
    lost = output_df[(output_df[sig_true]==1) & (output_df["Line_zero"] | output_df["Tester_zero"])]
    output_base = hse_file.split(".csv")[0] + "_putative_loss.csv"
    lost.to_csv(output_base, index=False)
    
    
if __name__ == "__main__":
    start = datetime.now()
    args = parse_arguments()
    stderr.write("Executed: python {0} -c {1} -b {2}\n".format(argv[0],args.COUNTS, args.BIAS))
    df = find_putative_loss(args.COUNTS)
    verify_HSE(df, args.BIAS)
    stop = datetime.now()
    stderr.write("Runtime: {0}\n".format(str(stop - start)))
