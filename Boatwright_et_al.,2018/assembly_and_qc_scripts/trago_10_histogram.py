#!/usr/bin/env python

from Bio import SeqIO
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


## Alireza Trinity for Tdu with 10 bp 
sizes = [len(rec) for rec in SeqIO.parse("/scratch/lfs/mcintyre/trago/outfiles/trinity_10/Tdu/Trinity.fasta", "fasta")]
len(sizes), min(sizes), max(sizes)
plt.hist(sizes, bins=10000)
plt.ylim([0,1600])
plt.xlim([0,5000])
plt.xlabel('seq length (only 0 to 5000)')
plt.ylabel('Tdu_10')
plt.savefig("/scratch/lfs/mcintyre/trago/outfiles/trinity_10/Tdu_10_plot")

## Alireza Trinity for Tpo with	10 bp 
sizes = [len(rec) for rec in SeqIO.parse("/scratch/lfs/mcintyre/trago/outfiles/trinity_10/Tpo/Trinity.fasta", "fasta")]
len(sizes), min(sizes), max(sizes)
plt.hist(sizes, bins=10000)
plt.ylim([0,1600])
plt.xlim([0,5000])
plt.xlabel('seq length (only 0 to 5000)')
plt.ylabel('Tpo_10')
plt.savefig("/scratch/lfs/mcintyre/trago/outfiles/trinity_10/Tpo_10_plot")

## Alireza Trinity for Tpr with	10 bp 
sizes = [len(rec) for rec in SeqIO.parse("/scratch/lfs/mcintyre/trago/outfiles/trinity_10/Tpr/Trinity.fasta", "fasta")]
len(sizes), min(sizes), max(sizes)
plt.hist(sizes,bins=10000)
plt.ylim([0,1600])
plt.xlim([0,5000])
plt.xlabel('seq length (only 0 to 5000)')
plt.ylabel('Tpr_10')
plt.savefig("/scratch/lfs/mcintyre/trago/outfiles/trinity_10/Tpr_10_plot")

