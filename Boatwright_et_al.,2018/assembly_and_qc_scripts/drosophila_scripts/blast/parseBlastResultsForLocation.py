#!/usr/bin/env python

# This script parses BLAST output to remove exact, self-hits.
# Output includes query, subject, percent identity, query start,
# query end, subject start, subject end, and the e-value in a 
# BED formatted file.
#AMR 03/28/2013

import csv
import operator

with open('/project/ambiguity/blast_results.tsv', 'rb') as input:
    with open('/project/ambiguity/blast_ambig_regions.bed', 'wb') as output:
        input_read = csv.reader(input, delimiter='\t')
        input_sort = sorted(input_read, key=operator.itemgetter(0)) # Sorts the input file by fusion_id
        for row in input_sort:
            query=row[0]
            subject=row[1]
            per_identity=row[2]
            q_start=int(row[6])-1 #BED files are 0-based, and BLAST results are 1-based, so 1 must be subtracted.
            q_end=row[7]
            s_start=int(row[8])+1 #BED files are 0-based, and BLAST results are 1-based, so 1 must be subtracted.
            s_end=row[9]
            e_value=row[10]
            if query==subject and per_identity=='100.00':
                continue
            else:
                output.write(query+'\t')
                #output.write(subject+'\t')
                #output.write(per_identity+'\t')
                output.write(str(q_start)+'\t')
                output.write(q_end+'\t')
                #output.write(str(s_start)+'\t')
                #output.write(s_end+'\t')
                output.write(e_value+'\n')


