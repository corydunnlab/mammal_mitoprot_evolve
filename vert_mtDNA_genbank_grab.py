#!/usr/bin/env python
# coding: utf-8 

# Funding received from the European Research Council, the Sigrid Jusélius Foundation, and the Academy of Finland contributed to the development of this software.
# Author: Cory Dunn
# Institution: University of Helsinki
# Author Email: cory.dunn@helsinki.fi
# License: GPLv3

# Load dependencies

import pandas as pd
import numpy as np
import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq

# Features set by user

genbank_file_to_use = 'Mammals_JAN_25_2021_and_A_punctatus.gb'
file_prefix = 'mammal_and_Ap_mtDNA_PYTEST_FEB_21_'
chosen_genes_set = {'ND1','ND2','COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5','ND6','CYTB'} # vertebrates

# Determine protein content and make lists

chosen_genes = list(chosen_genes_set)

CDS_list = []
gene_list = []
accession_error_list = []
accession_name_list = []
for seq_record in SeqIO.parse(genbank_file_to_use, "genbank"):
    accession = str(seq_record.id)
    accession_name_list.append(accession)
    for feat in range(len(seq_record.features)):
        feature_type = seq_record.features[feat].type
        if feature_type == 'CDS':
            try:
                prot_id = str(seq_record.features[feat].qualifiers['gene'][0])
                CDS_list.append(prot_id)
            except KeyError:
                accession_error_list.append(accession)
            
        if feature_type == 'gene':
            try:
                gene_id = str(seq_record.features[feat].qualifiers['gene'][0])
                gene_list.append(gene_id)
            except KeyError:
                accession_error_list.append(accession)
CDS_list_rem_dups = sorted(list(set(CDS_list)))
gene_list_rem_dups = sorted(list(set(gene_list)))
accession_error_list_rem_dups = sorted(list(set(accession_error_list)))
print('Protein set detected: ',CDS_list_rem_dups)

# Set up protein dataframe

protein_array = pd.DataFrame(index=accession_name_list,columns=CDS_list_rem_dups)
protein_array.index.name = 'Accession'

# Move translation products to protein dataframe

accession_name_list = []
for seq_record in SeqIO.parse(genbank_file_to_use, "genbank"):
    accession = str(seq_record.id)
    accession_name_list.append(accession)
    for feat in range(len(seq_record.features)):
        feature_type = seq_record.features[feat].type
        if feature_type == 'CDS':
            try:
                feature_translation = str(seq_record.features[feat].qualifiers['translation'][0])
                gene_id = str(seq_record.features[feat].qualifiers['gene'][0])
                protein_array.at[accession,gene_id] = feature_translation
            except KeyError:
                pass

# Count and record instances of each protein in the dataframe

protein_instances_list = []
for prot in CDS_list_rem_dups:
    #print(prot, len(protein_array)-protein_array[prot].isnull().sum())
    protein_to_append = prot, len(protein_array)-protein_array[prot].isnull().sum()
    protein_instances_list.append(protein_to_append) 

protein_instances_list_DF = pd.DataFrame(protein_instances_list,columns = ['Protein','Instances'])
protein_instances_list_DF.to_csv(file_prefix + 'protein_instances.csv',index=False)

#Set up gene dataframe

gene_array = pd.DataFrame(index=accession_name_list,columns=CDS_list_rem_dups)
gene_array.index.name = 'Accession'

# Move gene sequences to gene dataframe

accession_name_list = []
for seq_record in SeqIO.parse(genbank_file_to_use, "genbank"):
    accession = str(seq_record.id)
    accession_name_list.append(accession)
    gb_entire_sequence = list(seq_record.seq)
    gb_entire_sequence_lower = [x.lower() for x in gb_entire_sequence]
    gb_entire_sequence_joined = ''.join(gb_entire_sequence_lower)
    for feat in range(len(seq_record.features)):
        feature_type = seq_record.features[feat].type
        if feature_type == 'gene':
            try:
                feature_start_zero_based_numbering = seq_record.features[feat].location.nofuzzy_start
                feature_end_zero_based_numbering = seq_record.features[feat].location.nofuzzy_end
                feature_strand = seq_record.features[feat].strand
                gene_id = str(seq_record.features[feat].qualifiers['gene'][0])
                sequence_slice = gb_entire_sequence_joined[feature_start_zero_based_numbering:feature_end_zero_based_numbering]
                if feature_strand == 1:
                    gene_array.at[accession,gene_id] = sequence_slice
                if feature_strand == -1:
                    sequence_slice_BP = Seq(sequence_slice)
                    sequence_slice_rvscomp = str(sequence_slice_BP.reverse_complement())
                    gene_array.at[accession,gene_id] = sequence_slice_rvscomp
            except KeyError:
                pass

# Count and record instances of each gene in the dataframe

gene_instances_list = []
#print('Counts of each gene extracted from GenBank file:')
for gene in gene_list_rem_dups:
    gene_to_append = gene, len(gene_array)-gene_array[gene].isnull().sum()
    gene_instances_list.append(gene_to_append) 
    #print(gene, len(gene_array)-gene_array[gene].isnull().sum())
gene_instances_list_DF = pd.DataFrame(gene_instances_list,columns = ['Gene','Instances'])
gene_instances_list_DF.to_csv(file_prefix + 'gene_instances.csv',index=False)

# Generate and save a description dataframe

description_array = pd.DataFrame(index=accession_name_list,columns=['Organism'])
description_array.index.name = 'Accession'

accession_name_list = []
for seq_record in SeqIO.parse(genbank_file_to_use, "genbank"):
    accession = str(seq_record.id)
    for feat in range(len(seq_record.features)):
        feature_type = seq_record.features[feat].type
        if feature_type == 'source':
            try:
                organism = str(seq_record.features[feat].qualifiers['organism'][0])
                description_array.at[accession,'Organism'] = organism
            except KeyError:
                pass

description_array['Sequence_name'] = description_array.index + '_' + description_array['Organism']
description_array['Sequence_name'] = description_array['Sequence_name'].replace('[ .]', '_', regex=True)
description_array['Sequence_name'] = description_array['Sequence_name'].replace('[^a-zA-Z0-9_]', '', regex=True)

# Delete undesired columns based upon chosen genes

protein_array = protein_array[chosen_genes]
gene_array = gene_array[chosen_genes]

# Protein dataframe columns rename
    
protein_array_column_names_list = list(protein_array.columns.values.tolist())
protein_array_column_names_list_rename = [name + '_prot' for name in protein_array_column_names_list]
for i in range(len(protein_array_column_names_list)):
    protein_array.rename(columns = {protein_array_column_names_list[i]:protein_array_column_names_list_rename[i]}, inplace = True) 

# Dataframes concatenate sequences and form combined array

gene_array['Concatenate_Coding_DNA'] = ''

for gene in chosen_genes:
    gene_array['Concatenate_Coding_DNA'] = gene_array['Concatenate_Coding_DNA'] + gene_array[gene]

proteins_to_concatenate = []
for gene in chosen_genes:
    proteins_to_concatenate.append(gene + '_prot')

protein_array['Concatenate_Protein'] = ''
for protein in proteins_to_concatenate:
    protein_array['Concatenate_Protein'] = protein_array['Concatenate_Protein'] + protein_array[protein]
    
combined_array = description_array.join(protein_array).join(gene_array)

# Remove any entries from gene dataframe with empty cell

combined_array.replace('', np.nan, inplace=True)
combined_array.dropna(inplace=True)

# Remove duplicates from concatenates columns

combined_array.drop_duplicates(subset=['Concatenate_Coding_DNA'])
combined_array.drop_duplicates(subset=['Concatenate_Protein'])

# Choose sequence_name as index

combined_array.reset_index(drop=True, inplace=True)

# Save combined array

combined_array.to_csv(file_prefix + 'combined_input_array.csv', index=True)

print(str(len(combined_array)) + ' samples to be analyzed after pre-processing.')

# Saving gene concatenates

ofile = open(file_prefix+'coding_DNA_sequence_concatenate.fasta', "w")
for seqi in range(len(combined_array)):
    ofile.write(">" + combined_array.at[seqi,'Sequence_name'] + "\n" + combined_array.at[seqi,'Concatenate_Coding_DNA'] + "\n")
ofile.close()

# Saving protein concatenates and individual polypeptides

ofile = open(file_prefix+'protein_sequence_concatenate.fasta', "w")
for seqi in range(len(combined_array)):
    ofile.write(">" + combined_array.at[seqi,'Sequence_name'] + "\n" + combined_array.at[seqi,'Concatenate_Protein'] + "\n")
ofile.close()

for protcol in proteins_to_concatenate:
    ofile = open(file_prefix + protcol + '.fasta', "w")
    for seqi in range(len(combined_array)):
        ofile.write(">" + combined_array.at[seqi,'Sequence_name'] + "\n" + combined_array.at[seqi,protcol] + "\n")
    ofile.close()

# Free memory

del combined_array
del gene_array
del protein_array
