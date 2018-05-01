import pandas as pd
import sys

df = pd.read_table(sys.argv[1])
print(sys.argv[1])
print(len(df.index))
sig = df[df["adj.P.Val"] < 0.05]
not_sig = df[df["adj.P.Val"] >= 0.05]
print("Total significantly different (FDR 0.05): {0}".format(len(sig.index)))
print("Total not significantly different (FDR 0.05): {0}".format(len(not_sig.index)))

output_base = sys.argv[1].split('.')[0]
output_sig = output_base + "_sig_loci.txt"
output_not_sig = output_base + "_not_sig_loci.txt"

sig.to_csv(output_sig)
not_sig.to_csv(output_not_sig)
