# coding: utf-8

import pandas as pd
from glob import glob

glob("both_counts_T*")
tm = sorted([i for i in glob("both_counts_T*") if "Tdu_Tpo" in i])
tms = sorted([i for i in glob("both_counts_T*") if "Tdu_Tpr" in i])

col_name = "_".join(tm[0].split("_")[2:4])
df = pd.read_csv(tm[0])
df.rename(columns={"FUSION_ID":"FUSION_ID","BOTH_ALL":col_name}, inplace=True)
for file in tm[1:]:
    col_name = "_".join(file.split("_")[2:4])
    temp_df = pd.read_csv(file)
    temp_df.rename(columns={"FUSION_ID":"FUSION_ID","BOTH_ALL":col_name}, inplace=True)
    df = pd.merge(df, temp_df)

df.to_csv("both_count_Tdu_Tpo_all_reps.csv", index=False)



col_name = "_".join(tms[0].split("_")[2:4])
df = pd.read_csv(tms[0])
df.rename(columns={"FUSION_ID":"FUSION_ID","BOTH_ALL":col_name}, inplace=True)
for file in tms[1:]:
    col_name = "_".join(file.split("_")[2:4])
    temp_df = pd.read_csv(file)
    temp_df.rename(columns={"FUSION_ID":"FUSION_ID","BOTH_ALL":col_name}, inplace=True)
    df = pd.merge(df, temp_df)

df.to_csv("both_count_Tdu_Tpr_all_reps.csv", index=False)
