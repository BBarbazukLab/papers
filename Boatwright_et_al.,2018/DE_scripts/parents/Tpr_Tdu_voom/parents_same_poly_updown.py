# coding: utf-8
import pandas as pd

add_diff = pd.read_csv("Tms_sig_loci_parents_same.csv", index_col=0)
all_samples = pd.read_csv("all_samples.csv", index_col=0)

merged = add_diff.join(all_samples, how='right')

new = merged[~merged.t.isnull()]
new["Tms_mean"] = (new.Tms_1 + new.Tms_2 + new.Tms_3)/3
new["Tdu_mean"] = (new.Tdu_1 + new.Tdu_2 + new.Tdu_3)/3
new["Tpr_mean"] = (new.Tpr_1 + new.Tpr_2 + new.Tpr_3)/3

print(len(new.t))
print("Tms highest: {0}".format(len(new[(new.Tms_mean >= new.Tpr_mean) & (new.Tms_mean >= new.Tdu_mean)])))
print("Tms lowest: {0}".format(len(new[(new.Tms_mean <= new.Tpr_mean) & (new.Tms_mean <= new.Tdu_mean)])))
