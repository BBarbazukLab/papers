# coding: utf-8
import pandas as pd

add_diff = pd.read_csv("Tm_sig_loci_parents_same.csv", index_col=0)
all_samples = pd.read_csv("all_samples.csv", index_col=0)

merged = add_diff.join(all_samples, how='right')

new = merged[~merged.t.isnull()]
new["Tm_mean"] = (new.Tm_1 + new.Tm_2 + new.Tm_3)/3
new["Tdu_mean"] = (new.Tdu_1 + new.Tdu_2 + new.Tdu_3)/3
new["Tpo_mean"] = (new.Tpo_1 + new.Tpo_2 + new.Tpo_3)/3

print(len(new.t))
print("Tm highest: {0}".format(len(new[(new.Tm_mean >= new.Tpo_mean) & (new.Tm_mean >= new.Tdu_mean)])))
print("Tm lowest: {0}".format(len(new[(new.Tm_mean <= new.Tpo_mean) & (new.Tm_mean <= new.Tdu_mean)])))
