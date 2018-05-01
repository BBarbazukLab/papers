#!/usr/bin/env python

from Bio import SeqIO
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


## AMM normalized for Tdu no 10
sizes = [len(rec) for rec in SeqIO.parse("/scratch/lfs/mcintyre/trago/trago_output/trag_normalized_trinity_SC/Tdu_dir_denovTrinityOutput.Trinity.fasta", "fasta")]
len(sizes), min(sizes), max(sizes)
plt.hist(sizes, bins=10000)
plt.ylim([0,1600])
plt.xlim([0,5000])
plt.xlabel('seq length')
plt.ylabel('Tdu_amm_norm')
plt.savefig("/scratch/lfs/mcintyre/trago/trago_output/trag_normalized_trinity_SC/Tdu_amm_norm_plot2")

## AMM normalized for Tdu w/10
sizes = [len(rec) for rec in SeqIO.parse("/scratch/lfs/mcintyre/trago/trago_output/trag_normalized_trinity_SC/Tdu_10_dir_denovTrinityOutput.Trinity.fasta", "fasta")]
len(sizes), min(sizes), max(sizes)
plt.hist(sizes, bins=10000)
plt.ylim([0,1600])
plt.xlim([0,5000])
plt.xlabel('seq length')
plt.ylabel('Tdu_10_amm_norm')
plt.savefig("/scratch/lfs/mcintyre/trago/trago_output/trag_normalized_trinity_SC/Tdu_10_amm_norm_plot2")


## AMM normalized for Tpo no 10                    
sizes = [len(rec) for rec in SeqIO.parse("/scratch/lfs/mcintyre/trago/trago_output/trag_normalized_trinity_SC/Tpo_dir_denovTrinityOutput.Trinity.fasta", "fasta")]
len(sizes), min(sizes), max(sizes)
plt.hist(sizes, bins=10000)
plt.ylim([0,1600])
plt.xlim([0,5000])
plt.xlabel('seq length')
plt.ylabel('Tpo_amm_norm')
plt.savefig("/scratch/lfs/mcintyre/trago/trago_output/trag_normalized_trinity_SC/Tpo_amm_norm_plot2")

## AMM normalized for Tpo w/10
sizes = [len(rec) for rec in SeqIO.parse("/scratch/lfs/mcintyre/trago/trago_output/trag_normalized_trinity_SC/Tpo_10_dir_denovTrinityOutput.Trinity.fasta", "fasta")]
len(sizes), min(sizes), max(sizes)
plt.hist(sizes, bins=10000)
plt.ylim([0,1600])
plt.xlim([0,5000])
plt.xlabel('seq length')
plt.ylabel('Tpo_10_amm_norm')
plt.savefig("/scratch/lfs/mcintyre/trago/trago_output/trag_normalized_trinity_SC/Tpo_10_amm_norm_plot2")

## AMM normalized for Tpr no 10
sizes = [len(rec) for rec in SeqIO.parse("/scratch/lfs/mcintyre/trago/trago_output/trag_normalized_trinity_SC/Tpr_dir_denovTrinityOutput.Trinity.fasta", "fasta")]
len(sizes), min(sizes), max(sizes)
plt.hist(sizes, bins=10000)
plt.ylim([0,1600])
plt.xlim([0,5000])
plt.xlabel('seq length')
plt.ylabel('Tpr_norm')
plt.savefig("/scratch/lfs/mcintyre/trago/trago_output/trag_normalized_trinity_SC/Tpr_amm_norm_plot2")

## AMM normalized for Tpr w/10
sizes = [len(rec) for rec in SeqIO.parse("/scratch/lfs/mcintyre/trago/trago_output/trag_normalized_trinity_SC/Tpr_10_dir_denovTrinityOutput.Trinity.fasta", "fasta")]
len(sizes), min(sizes), max(sizes)
plt.hist(sizes, bins=10000)
plt.ylim([0,1600])
plt.xlim([0,5000])
plt.xlabel('seq length')
plt.ylabel('Tpr_10_norm')
plt.savefig("/scratch/lfs/mcintyre/trago/trago_output/trag_normalized_trinity_SC/Tpr_10_amm_norm_plot2")

