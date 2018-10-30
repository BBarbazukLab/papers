#!/bin/bash

# bash generate_regulator_networks.bash regulatory_effects_0-7.txt
grep -P "^\d+" $1 | cut -f 5,7 | tr ',' '\n' | tr '\t' '\n' | sort -u > regulatory_genes.txt

python robust_filter_coexpression_gene_pairs.py regulatory_genes.txt

if [ -e regulatory_genes.counts.txt ]
  then
    rm regulatory_genes.counts.txt
fi

while read p
  do
    echo "${p}" >> regulatory_genes.counts.txt
    cut -f 6 -d ',' regulatory_genes.csv | grep -ciP "^${p}_|_${p}$" >> regulatory_genes.counts.txt
  done < regulatory_genes.txt

grep -B1 -P "^[3-9]|^[1-9][0-9]+" regulatory_genes.counts.txt | grep -P "^\D+" | grep -v "\-\-" > ge_three_edges.txt
grep -B1 -P "^[5-9]|^[1-9][0-9]+" regulatory_genes.counts.txt | grep -P "^\D+" | grep -v "\-\-" > ge_five_edges.txt
grep -B1 -P "^[1-9][0-9]+" regulatory_genes.counts.txt | grep -P "^\D+" | grep -v "\-\-" > ge_ten_edges.txt

while read p
  do
    # Remove csv if exists
    if [ -e ${p}.csv ]
      then
        rm ${p}.csv
    fi

    # For each gene, generate an interaction cluster
    echo "SUID,interaction,name,selected,shared interaction,shared name,SP,A,B" > ${p}.csv
    grep -i "${p}" regulatory_genes.csv >> ${p}.csv

    # For each gene.csv cluster -- generate table of nodes
    python genes_from_edge_table.py ${p}.csv ${p} > ${p}.txt

    # generate plot for table of nodes
    python pointplot.py aco_vs_mus_counts.txt ${p}.txt

    # Clean up
    mv figure.pdf ${p}.pdf
#    rm ${p}.csv 
#    rm ${p}.txt

  done < ge_ten_edges.txt
