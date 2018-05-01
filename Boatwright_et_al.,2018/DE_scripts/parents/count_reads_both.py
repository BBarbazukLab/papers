import pandas as pd
import sys

if len(sys.argv) != 2:
    print("python {0} ase_counts_Tdu_1_2_Tdu_Tpo.csv".format(sys.argv[0]))
    sys.exit(1)

in_name  = sys.argv[1]
out_name = "both_" + "_".join(in_name.split("_")[1:])

#              0        1            2                    3               4                     5                          6                        7                           8                    9                   10
# HEADER: FUSION_ID,BOTH_EXACT,BOTH_INEXACT_EQUAL,SAM_A_ONLY_EXACT,SAM_B_ONLY_EXACT,SAM_A_EXACT_SAM_B_INEXACT,SAM_B_EXACT_SAM_A_INEXACT,SAM_A_ONLY_SINGLE_INEXACT,SAM_B_ONLY_SINGLE_INEXACT,SAM_A_INEXACT_BETTER,SAM_B_INEXACT_BETTER
df = pd.read_csv(in_name)
df["BOTH_ALL"]=df[["BOTH_EXACT","BOTH_INEXACT_EQUAL","SAM_A_EXACT_SAM_B_INEXACT","SAM_B_EXACT_SAM_A_INEXACT","SAM_A_INEXACT_BETTER","SAM_B_INEXACT_BETTER"]].sum(axis=1)

df.rename(columns={"BOTH_ALL":"Count"}, inplace=True)
df[["FUSION_ID","Count"]].to_csv(out_name, index=False)
