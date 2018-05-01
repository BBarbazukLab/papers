import pandas as pd
from collections import OrderedDict

# Last file is the additive one of interest
files = ["Tm_sig_loci_parents_same.csv","Tm_sig_loci_parents_diff.csv","Tm_not_sig_parents_same.csv","Tm_not_sig_parents_diff.csv"]

loci = []
factor_labels = OrderedDict()
for file in files:
    df = pd.read_csv(file)
    loci += list(df["Unnamed: 0"])
    if file == "Tm_not_sig_parents_diff.csv":
        for locus in df["Unnamed: 0"]:
            factor_labels[locus]="Tm_additive"

with open("enrichment_loci.txt",'w') as f:
    for i in loci:
        f.write(i + "\n")

with open("factor_labeling.txt",'w') as f:
    for key, value in factor_labels.items():
        f.write(value + "\t" + key + "\n")
