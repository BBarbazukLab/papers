#!/usr/bin/env python

from Bio import SeqIO
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


## Alireza Trinity for Tdu
sizes = [len(rec) for rec in SeqIO.parse("/scratch/lfs/mcintyre/trago/outfiles/trinity/Tdu/Trinity.fasta", "fasta")]
len(sizes), min(sizes), max(sizes)
plt.hist(sizes, bins=10000)
plt.ylim([0,1600])
plt.xlim([0,5000])
plt.xlabel('seq length (only 0 to 5000)')
plt.ylabel('Tdu')
plt.savefig("/scratch/lfs/mcintyre/trago/outfiles/trinity/Tdu_plot")

## Alireza Trinity for Tpo
sizes = [len(rec) for rec in SeqIO.parse("/scratch/lfs/mcintyre/trago/outfiles/trinity/Tpo/Trinity.fasta", "fasta")]
len(sizes), min(sizes), max(sizes)
plt.hist(sizes, bins=10000)
plt.ylim([0,1600])
plt.xlim([0,5000])
plt.xlabel('seq length (only 0 to 5000)')
plt.ylabel('Tpo')
plt.savefig("/scratch/lfs/mcintyre/trago/outfiles/trinity/Tpo_plot")

## Alireza Trinity for Tpr
sizes = [len(rec) for rec in SeqIO.parse("/scratch/lfs/mcintyre/trago/outfiles/trinity/Tpr/Trinity.fasta", "fasta")]
len(sizes), min(sizes), max(sizes)
plt.hist(sizes,bins=10000)
plt.ylim([0,1600])
plt.xlim([0,5000])
plt.xlabel('seq length (only 0 to 5000)')
plt.ylabel('Tpr')
plt.savefig("/scratch/lfs/mcintyre/trago/outfiles/trinity/Tpr_plot")

