# Rename Trinity.fasta -- identifiers are not compatible with WU-BLAST...
python ~/scripts/rename_FASTA_deflines.py -f Trinity.fasta -d Trinity_contig
# generates:
Trinity.renamed.table
Trinity.renamed.fa

# Use WU-BLAST to Mus-to-Acomys and Acomys-to-Mus

# Parse the blastn output
# no parameters set other than nextQuery or not (XY_WU_BP_mod_next.pl or XY_WU_BP_mod_no_next.pl)

# WU-BLAST output example:
parsed_Trinity_to_Mus.blastn.all
parsed_Mus_to_Trinity.blastn.all

# Find reciprocated hits
# executed from: PARSED_BLAST
python get_reciprocated.py -a parsed_Mus_to_Trinity.blastn.all -b parsed_Trinity_to_Mus.blastn.all
# generates:
reciprocated_blast_hits.txt
reciprocated_blast_hits.gene_trans.txt

# Filter reciprocated_blast_hits
python filter_reciprocated_blast.py -m reciprocated_blast_hits.gene_trans.txt -a parsed_Mus_to_Trinity.blastn.all 2> not_reciprocated.txt
# generates:
parsed_Mus_to_Trinity.blastn.all.filtered.txt

# Regenerate reciprocated hits using quality filtered reciprocal hits
cut -f 1,3 parsed_Mus_to_Trinity.blastn.all.filtered.txt | sort -u > quality_reciprocated_blast_hits.gene_trans.txt
# generates:q
quality_reciprocated_blast_hits.gene_trans.txt

# Generate BED files
python generate_BED_from_Mus_BLAST_gene_trans_map.py -m quality_reciprocated_blast_hits.gene_trans.txt -a parsed_Mus_to_Trinity.blastn.all.filtered.txt
# generates:
parsed_Mus_to_Trinity.blastn.all.filtered_Q.bed
parsed_Mus_to_Trinity.blastn.all.filtered_S.bed

# Correct Trinity names back to the original
python correct_renamed_contigs_in_BED.py Trinity_transmap.txt Trinity.renamed.table parsed_Mus_to_Trinity.blastn.all.filtered_S.renamed.bed
# generates:
parsed_Mus_to_Trinity.blastn.all.filtered_S.corrected.bed

# Merge overlapping coordinates in BED files
module add bedtools 
bedtools sort -i parsed_Mus_to_Trinity.blastn.all.filtered_Q|S.corrected.bed > parsed_Mus_to_Trinity.blastn.all.filtered_Q|S.sorted.bed
bedtools merge -i parsed_Mus_to_Trinity.blastn.all.filtered_Q|S.sorted.bed > parsed_Mus_to_Trinity.blastn.all.filtered_Q|S.sorted.merged.bed
# generates:
parsed_Mus_to_Trinity.blastn.all.filtered_Q|S.sorted.merged.bed

# Correct contig names for Trinity contigs in Transmap
python correct_renamed_contigs_2.py quality_reciprocated_blast_hits.gene_trans.txt Trinity.renamed.table 
# generates:
quality_reciprocated_blast_hits_corrected.txt

# Change Mus rna names to gene symbols
python correct_renamed_contigs_Mus_rna_to_gene.py quality_reciprocated_blast_hits_corrected.txt name_conversion.table 
# generates:
quality_reciprocated_blast_hits_corrected_corrected.txt --> renamed to Acomys-Mus_gene_transmap.txt

# Generate full transmap -- this script takes a few minutes to complete
python generate_full_transmap.py Trinity_transmap.txt quality_reciprocated_blast_hits_corrected_corrected.txt > Acomys_non-orthologous_gene_transmap.txt
# generates:
Acomys_non-orthologous_gene_transmap.txt

# Expand original gene trans with recovered isoforms
cat Acomys-Mus_gene_transmap.txt Acomys_non-orthologous_gene_transmap.txt | sort > Acomys_per_se_gene_transmap.txt
# generates:
Acomys_per_se_gene_transmap.txt

# Correct Trinity contig and Mus rna names in BLAST (optional for visualization)
python correct_renamed_contigs_in_BLAST.py parsed_Mus_to_Trinity.blastn.all.filtered.txt name_conversion.table Trinity.renamed.table
# generates:
parsed_Mus_to_Trinity_corrected.txt
