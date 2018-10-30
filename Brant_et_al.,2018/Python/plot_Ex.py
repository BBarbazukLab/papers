# coding: utf-8

"""plot_Ex.py: This script generates a scatter plot from the data generated by Trinity's contig_ExN50_statistic.pl script. This script will write the plot to Ex_vs_ExN50.png by default."""

__author__ = "J. Lucas Boatwright"

import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("ExN50.stats",sep="\t")
df['E']=df['#E'].apply(lambda x: int(x.lstrip('E')))

df.plot.scatter('E','E-N50')
plt.xlabel("Ex")
plt.ylabel("Ex-N50")
plt.savefig("Ex_vs_ExN50")
