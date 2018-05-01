import pandas as pd
import numpy as np
from sys import argv, exit

if len(argv) != 2:
    print("python {0} ase_bayes_flagged_counts.csv".format(argv[0]))
    exit(1)

def fix_flags():
    """Fix the flags so that sum > 0 == 1 else 0"""
    df = pd.read_csv(argv[1], index_col=0)
    df.pop("flag_analyze")
    df["flag_analyze"]=np.where(df.sum(axis=1)>0,1,0)
    output = argv[1].split(".csv")[0] + "_re-flagged.csv"
    df.to_csv(output)


if __name__ == "__main__":
    fix_flags()
