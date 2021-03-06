# Documentation

Read pre-processing/QC/assembly construction are described in worksheet_trago_amm5.xls

Read alignment/Poisson-Gamma/final result processing may be viewed in <a href="https://htmlpreview.github.io/?https://github.com/BBarbazukLab/papers/blob/master/Boatwright_et_al.%2C2018/New_World_PG_pipeline_documentation.html" target="_blank">New_World_PG_pipeline_documentation.html</a>.

# Python_CORE_scripts
All python scripts contain usage information that may be accessed using:

```bash
python <script_name> -h
```

The help usage is as follows:

```bash
python group_self_BLAST.py -h

usage: group_self_BLAST.py [-h] -b BLAST

This script is designed to take a parsed BLAST file and generate a     
table of groups. The grouping method is greedy.

     Example: python group_self_BLAST.py -b parsed_BLAST

optional arguments:
  -h, --help            show this help message and exit

required arguments:
  -b BLAST, --BLAST BLAST
                        Parsed BLAST file used for grouping transcripts.
```

```bash
python contig_compression.py -h

usage: contig_compression.py [-h] -f FASTA

This script takes each group from grouped_hits.txt and executes     
CAP3 on each group in an attempt to merge transcripts and remove     
redundant contigs. Requires a FASTA file containing all sequences from     
grouped_hits.txt. Auto-detects 'grouped_hits.txt' in the current     
directory. All contigs after compression are in 'contigs_collapsed.fasta'

     Example: python contig_compression.py -d contigs.fa

optional arguments:
  -h, --help            show this help message and exit

required arguments:
  -f FASTA, --FASTA FASTA
                        FASTA file containing all sequences within
                        grouped_hits.txt.
```

```bash
python isolate_consensus.py -h

usage: isolate_consensus.py [-h] -g GROUPS -f FASTA

This script receives 'grouped_hits.txt', uses a FASTA file to get the     
sequences, uses MAFFT to align each group, then generates a consensus     
sequence for each group using biopython.

     Example: python isolate_consensus.py -g grouped_hits.txt -f sequences.fasta

optional arguments:
  -h, --help            show this help message and exit

required arguments:
  -g GROUPS, --GROUPS GROUPS
                        Grouped hits file generated by group_self_BLAST.py
  -f FASTA, --FASTA FASTA
                        FASTA file containing all sequences.
```

```bash
python BLAST_COREs.py -h

usage: BLAST_COREs.py [-h] -a BLAST_A -b BLAST_B -p PROGRAM

This script is designed to identify reciprocal best-hits given two     
parsed BLAST files. Output includes reciprocated_blast_hits.txt and     
two BED files representing HSP coordinates from BLAST_A for query     
and subject.
    
    Usage: python BLAST_COREs.py -a A_to_B_BLAST -b B_to_A_BLAST    

optional arguments:
  -h, --help            show this help message and exit
  -a BLAST_A, --BLAST_A BLAST_A
                        Parsed BLAST file from species A in tab-separated
                        format.
  -b BLAST_B, --BLAST_B BLAST_B
                        Parsed BLAST file from species B in tab-separated
                        format.
  -p PROGRAM, --PROGRAM PROGRAM
                        Program used to generate tabular blast -- NCBI or WU
                        parsed.
```

```bash
python determine_COREs.py -h

usage: determine_COREs.py [-h] {known,rbh} ...

This script is designed to identify Common Orthologous REgions (COREs)     
given either two parsed and filtered BLAST files and the FASTA files     
containing the corresponding sequences (see 'rbh' subcommand) or given a     
CSV with orthologous pairs and the FASTA files containing the     
corresponding sequences (see 'known' subcommand). If 'rbh' is used,     
BLAST files are assumed to have a single best hit for each query.     
COREs are found in Species<AB>-species<BA>_ortholog_COREs.bed.
    
    Example (known): python determine_COREs.py -o orthologs.csv -f A.fasta -g B.fasta    
    Example (rbh): python determine_COREs.py -a parsed_BLAST_A -b parsed_BLAST_B -f A.fasta -g B.fasta    
    
Example BLAST_A:     
    contig_10340    101674  contig_252388   66445   1621    0.  80.4    80.4878048780488    533 32644   32124   65919   66445   -    

Example BLAST_B:     
    contig_252388   66865   contig_10340    101674  280762  0.  91.6    91.6699410609037    66170   66192   58  32362   98486   -    

positional arguments:
  {known,rbh}  help for subcommand
    known      If orthologous pairs are known, use a user-defined CSV with
               Species_A_GeneID,Species_B_GeneID
    rbh        If orthologs are not yet identified, use a reciprocal best-hit
               approach to determine orthologous pairs.

optional arguments:
  -h, --help   show this help message and exit
```

```bash
python3 BED_mean_difference.py -h

usage: BED_mean_difference.py [-h] -a BED_A -b BED_B [-o OUTPUT]

This script is designed to take two bed files and generate a mean-difference plot of the CORE lengths

     Example: python BED_mean_difference.py -a A.bed -b B.bed -o BED_mean_difference.pdf

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --OUTPUT OUTPUT
                        Output file name.

required arguments:
  -a BED_A, --BED_A BED_A
                        BED file A
  -b BED_B, --BED_B BED_B
                        BED file B
```

```bash
python3 BED_GC_MD.py -h

usage: BED_GC_MD.py [-h] -a BED_A_STATS -b BED_B_STATS -o ORTHOLOGS

This script is designed to generate an MD plot of GC bias

     Example: python BED_GC_MD.py -a A.bed.stats.csv -b B.bed.stats.csv -o orthlogs.txt

optional arguments:
  -h, --help            show this help message and exit

required arguments:
  -a BED_A_STATS, --BED_A_STATS BED_A_STATS
                        BED file A stats from pyfasta info
  -b BED_B_STATS, --BED_B_STATS BED_B_STATS
                        BED file B stats from pyfasta info
  -o ORTHOLOGS, --ORTHOLOGS ORTHOLOGS
                        Orthologs.
```
