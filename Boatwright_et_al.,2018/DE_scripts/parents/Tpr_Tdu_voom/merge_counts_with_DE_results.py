# coding: utf-8
import pandas as pd

add_diff = pd.read_csv("Tms_not_sig_parents_diff.csv", index_col=0)
all_samples = pd.read_csv("all_samples.csv",index_col=0)

merged = add_diff.join(all_samples, how='right')

new = merged[~merged.t.isnull()]
new["Tdu_mean"] = (new.Tdu_1 + new.Tdu_2 + new.Tdu_3)/3
new["Tpr_mean"] = (new.Tpr_1 + new.Tpr_2 + new.Tpr_3)/3

print(len(new.t))
print("Tdu higher: {0}".format(len(new[new.Tdu_mean > new.Tpr_mean])))
print("Tpr higher: {0}".format(len(new[new.Tdu_mean < new.Tpr_mean])))
