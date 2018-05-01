# coding: utf-8
import pandas as pd
df = pd.read_table("DE_Tdu_Tpo.txt")
same = pd.read_csv("DE_Tdu_Tpo_nonsigFDR.csv")
same_set = set(same["Unnamed: 0"])
df_set = set(df.index)
different = df_set.difference(same_set)
output = list(different)
with open("parents_different.txt",'w') as f:
    for line in output:
        f.write(line + "\n")
        
