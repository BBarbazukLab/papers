#!/usr/bin/env python

from Bio import SeqIO
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


## AMM uniq for Tdu no 10
sizes = [len(rec) for rec in SeqIO.parse("/scratch/lfs/mcintyre/trago/trago_output/trinity/trinity/Tdu.Trinity.fasta", "fasta")]
len(sizes), min(sizes), max(sizes)
plt.hist(sizes, bins=10000)
plt.ylim([0,1600])
plt.xlim([0,5000])
plt.xlabel('seq length')
plt.ylabel('Tdu_uniq_norm')
plt.savefig("/scratch/lfs/mcintyre/trago/trago_output/trinity/trinity/Tdu_amm_uniq_plot2")

## AMM uniq for Tdu no 10
sizes = [len(rec) for rec in SeqIO.parse("/scratch/lfs/mcintyre/trago/trago_output/trinity/trinity/Tdu.Trinity.fasta", "fasta")]
len(sizes), min(sizes), max(sizes)
plt.hist(sizes, bins=10000)
##plt.ylim([0,1600])
##plt.xlim([0,5000])
plt.xlabel('seq length')
plt.ylabel('Tdu_uniq_norm')
plt.savefig("/scratch/lfs/mcintyre/trago/trago_output/trinity/trinity/Tdu_amm_uniq_plot")



## AMM uniq for Tpo no 10
sizes = [len(rec) for rec in SeqIO.parse("/scratch/lfs/mcintyre/trago/trago_output/trinity/trinity/Tpo.Trinity.fasta", "fasta")]
len(sizes), min(sizes), max(sizes)
plt.hist(sizes, bins=10000)
plt.ylim([0,1600])
plt.xlim([0,5000])
plt.xlabel('seq length')
plt.ylabel('Tpo_uniq_norm')
plt.savefig("/scratch/lfs/mcintyre/trago/trago_output/trinity/trinity/Tpo_amm_uniq_plot2")

## AMM uniq for Tpo no 10
sizes = [len(rec) for rec in SeqIO.parse("/scratch/lfs/mcintyre/trago/trago_output/trinity/trinity/Tpo.Trinity.fasta", "fasta")]
len(sizes), min(sizes), max(sizes)
plt.hist(sizes, bins=10000)
##plt.ylim([0,1600])
##plt.xlim([0,5000])
plt.xlabel('seq length')
plt.ylabel('Tpo_uniq_norm')
plt.savefig("/scratch/lfs/mcintyre/trago/trago_output/trinity/trinity/Tpo_amm_uniq_plot")



## AMM uniq for Tpr no 10
sizes = [len(rec) for rec in SeqIO.parse("/scratch/lfs/mcintyre/trago/trago_output/trinity/trinity/Tpr.Trinity.fasta", "fasta")]
len(sizes), min(sizes), max(sizes)
plt.hist(sizes, bins=10000)
plt.ylim([0,1600])
plt.xlim([0,5000])
plt.xlabel('seq length')
plt.ylabel('Tpr_uniq_norm')
plt.savefig("/scratch/lfs/mcintyre/trago/trago_output/trinity/trinity/Tpr_amm_uniq_plot2")

## AMM uniq for Tpr no 10
sizes = [len(rec) for rec in SeqIO.parse("/scratch/lfs/mcintyre/trago/trago_output/trinity/trinity/Tpr.Trinity.fasta", "fasta")]
len(sizes), min(sizes), max(sizes)
plt.hist(sizes, bins=10000)
##plt.ylim([0,1600])
##plt.xlim([0,5000])
plt.xlabel('seq length')
plt.ylabel('Tpr_uniq_norm')
plt.savefig("/scratch/lfs/mcintyre/trago/trago_output/trinity/trinity/Tpr_amm_uniq_plot")

