
"""generate_factor_plot_df.py: This script generates a dataframe by extracting genes from the expression matrix. This dataframe is then intended for use downstream by factor_plot.py"""

__author__      = "J. Lucas Boatwright"

import pandas as pd
import numpy as np

def generate_df(df, genes):
    "Select genes from df, rotate for plotting"
    # ex. genes = ["0610007P14Rik","0610009O20Rik","0610012G03Rik"]
    plotting_df = pd.DataFrame()
    df.index = pd.Series(df.index).apply(lambda x: x.lower())
    for gene in genes:
        gene_df = df.ix[gene]
        gene_df = gene_df.reset_index()
        gene_df = gene_df.drop(0)
        series = pd.Series([gene for i in range(len(df.columns)-1)], name="Gene Symbol", index=gene_df.index)        
        gene_df["index"] = gene_df["index"].apply(lambda x: x.split('.')[0])
        gene_df["Gene Symbol"] = series
        gene_df.columns = ["Sample","Expression","Gene Symbol"]
        gene_df["Expression"] = gene_df.Expression + 1.0
        gene_df["log2(Expression+1)"] = np.log2(gene_df["Expression"].astype(np.float))
        try:
            plotting_df = pd.concat([plotting_df, gene_df])
        except pd.tools.merge.MergeError:
            plotting_df = gene_df 
    return plotting_df
