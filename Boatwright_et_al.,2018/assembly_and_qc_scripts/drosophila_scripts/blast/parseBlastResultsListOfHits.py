#!/usr/bin/env python

# This script parses out a BLAST table to give out 
# query, number of BLAST hits, and a list of the BLAST hits. 
# Here, the script is pulling out dmel fusions that are
# similar to other dmel fusions. 

# Revision: The script now includes instances when 
# BLAST identifies a fusion hitting within itself
# or hitting with another fusion in multiple places.

# Created 02/20/2012 by Alexi Runnels and Chelsea Tymms

input = '/project/ambiguity/blast_results.tsv'
output= '/project/ambiguity/parsed_blast_list_of_hits.csv'
with open (input, 'rb') as file:
    a=0
    array=[] 
    newarray=[]
    for line in file:
        array.append([])
#split the BLAST file by tabs into first column (query), 
#second column(hits), and everything else
        query, hits, other = line.split('\t', 2)
        array[a].append(query)
        array[a].append(hits)
        a=a+1
    
    #now your array is filled
    i=0
    while i<len(array):
        str = array[i][0]
        j=0
        count=0
        newstr=[]
        while (i+j)<len(array) and array[i+j][0]==str:
            if str != array[i+j][1]:
                newstr.append(array[i+j][1])
            else:
                count += 1
                #print count
            j=j+1
        if count > 1:
            newstr.append(str)
        #newstr_no_dups=sorted(set(newstr)) #remove duplicates #This line is necessary if want only non-self hits. It is excluded here so that self hits other than those complete, exact self-matches are included in the results. 
        size=len(newstr)
        #for x in newstr_no_dups:
        #    x.replace('\'', '')
        #    print x
        #print only those matching non-self hits
        if size !=0:
            newarray.append([str, size, newstr])
        i=i+j
outfile = open(output, 'wb')
for entry in newarray:
    print >>outfile, entry
