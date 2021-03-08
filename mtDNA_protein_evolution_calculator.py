#!/usr/bin/env python
# coding: utf-8 

# Funding received from the European Research Council, the Sigrid Jusélius Foundation, and the Academy of Finland contributed to the development of this software.
# Authors: Cory Dunn and Bala Ani Akpinar
# Institution: University of Helsinki
# Author Email: cory.dunn@helsinki.fi
# License: GPLv3

# Load dependencies

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import os
import sys

from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
from sklearn.preprocessing import MinMaxScaler
from sklearn.utils import resample
from sklearn.metrics import accuracy_score

from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo
from Bio import Entrez

import csv
from itertools import zip_longest

# Features set by user

Entrez.api_key = "<INCLUDE USER ENTREZ KEY HERE>"
Entrez.email = "<INCLUDE USER EMAIL>"

build_tree_on = 'provided' # blank on provided tree, otherwise 'aminoacid' or 'nucleotide'
provided_tree = '/Users/corydunn/Dropbox/Lab/Current_Lab_Members/CDUNN/University_of_Helsinki/CDD_Mammal_Selection_JAN_25_2021/mammal_tree_building_JAN_25_2021/T3_MAMM_JAN_25_raxml_bestTree_BLonly_rooted_Ap.nwk'
taxonomy_groups_to_test = ['Class','Order','Suborder','Infraorder','Family','Subfamily','Genus']

genbank_file_to_use = '/Users/corydunn/Dropbox/Lab/Current_Lab_Members/CDUNN/University_of_Helsinki/CDD_Mammal_Selection_JAN_25_2021/Mammals_JAN_25_2021_and_A_punctatus.gb'
selected_accession = 'NC_006853_1_Bos_taurus'
file_prefix = 'mammal_and_Ap_mtDNA_PYTEST_FEB_21_'
chosen_genes_set = {'ND1','ND2','COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5','ND6','CYTB'} # vertebrates

gap_cutoff = 0.02

arctan_columns_to_assess = ['Sum-of-protein-scores_arctan','Complex-I_arctan','Complex-III_arctan','Complex-IV_arctan','Complex-V_arctan'] # Output columns of interest
alpha = 0.90 # Confidence level

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

print(str(len(combined_array)) + ' samples to be analyzed after pre-processing.')

# Free memory

del combined_array
del gene_array
del protein_array

# Run MAFFT with FFT_NS_2 on the output FASTA files

for protein_choice in proteins_to_concatenate:
    os.system("mafft --thread 4 --retree 2 --inputorder " + file_prefix + protein_choice + ".fasta > " + file_prefix + protein_choice + "_MAFFT_FFT_NS_2.fasta")

# PAGAN ancestral positions based upon NWK tree generated above

if build_tree_on == 'nucleotide':
    newick_tree = 'coding_DNA_sequence_concatenate_MAFFT_FFT_NS_2.nwk'
elif build_tree_on == 'aminoacid':
    newick_tree = 'protein_sequence_concatenate_MAFFT_FFT_NS_2.nwk'
elif build_tree_on == 'provided':
    newick_tree = provided_tree
    
for protein in chosen_genes:
    os.system("pagan -s " + file_prefix + protein + '_prot_MAFFT_FFT_NS_2.fasta' + " -t " + newick_tree + " -o " + file_prefix + protein + \
             "_prot_PAGAN --output-ancestors --guidetree --threads 6")

# Ungap PAGAN alignments based on human reference sequence

selected_PAGAN_alignments = []

for protein in chosen_genes:
    appendage =  file_prefix + protein + '_prot_PAGAN.fas'
    selected_PAGAN_alignments.append(appendage)

for alignfile in selected_PAGAN_alignments:

    record_x_toward_seq_dataframe = []
    sequence_records = []
    alignment_record_name_list = []

    for record in SeqIO.parse(alignfile,"fasta"):
        alignment_record_name_list.append(record.name)
        record_x_toward_seq_dataframe = list(record.seq)
        record_x_toward_seq_dataframe_UPPER = [x.upper() for x in record_x_toward_seq_dataframe] 
        sequence_records.append(record_x_toward_seq_dataframe_UPPER)

    depth_of_alignment = (len(alignment_record_name_list))

    accession_name_dataframe = pd.DataFrame(alignment_record_name_list, columns=['Accession'])
    sequence_dataframe = pd.DataFrame(sequence_records)
    
    
    my_accession_row = accession_name_dataframe[accession_name_dataframe['Accession'] == selected_accession].index[0]
    my_accession_sequence = sequence_dataframe.iloc[my_accession_row,:]

    sequence_dataframe_series_remove = my_accession_sequence == '-'
    sequence_dataframe_index_remove = sequence_dataframe_series_remove[sequence_dataframe_series_remove].index
    sequence_dataframe = sequence_dataframe.drop(sequence_dataframe_index_remove,axis=1)
    
    sequence_dataframe_concat = sequence_dataframe.apply(''.join, axis=1)
    sequence_dataframe_final_seq_list = sequence_dataframe_concat.values.tolist()
    sep = '.'
    filestem = alignfile.split(sep, 1)[0]
    ofile = open(filestem + '_' + selected_accession + '.fasta', "w")
    for seqi in range(len(alignment_record_name_list)):
        ofile.write(">" + alignment_record_name_list[seqi] + "\n" + sequence_dataframe_final_seq_list[seqi] + "\n")
    ofile.close()

# Extract tree distances and prepare a dataframe with edge names

newick_tree_to_use = file_prefix + chosen_genes[0]+'_prot_PAGAN.anctree'
distances_as_BL = []

def extract_distances(tree):

    edge_distances = {}

    for parent in tree.find_clades(terminal=False, order='level'):
        for child in parent.clades:
            cum_dist = 0
            for stop in tree.get_path(child):
                cum_dist += stop.branch_length

            edge_distances[parent.name + '*' + child.name] = child.branch_length,cum_dist

    return edge_distances

# Read the input tree

my_tree = Phylo.read(newick_tree_to_use, 'newick')

distances_as_BL = extract_distances(my_tree)

recip_TSS_plus_one_proteins_SUMS_BL = pd.DataFrame(distances_as_BL).transpose()
recip_TSS_plus_one_proteins_SUMS_BL.columns =['Branch_length','Distance_from_root']
recip_TSS_plus_one_proteins_SUMS_BL.index.name = 'Edge_name'
recip_TSS_plus_one_proteins_SUMS_BL.reset_index(inplace=True)

# Begin taxonomy recovery

taxonomy_array = description_array.copy()
taxonomy_columns_A = ['TaxonID'] + taxonomy_groups_to_test
taxonomy_columns_B = taxonomy_groups_to_test
taxonomy_array[taxonomy_columns_A] = ""

# Move taxon information sequences to taxon dataframe

for seq_record in SeqIO.parse(genbank_file_to_use, "genbank"):
    accession = str(seq_record.id)
    for feat in range(len(seq_record.features)):
        feature_type = seq_record.features[feat].type
        if feature_type == 'source':
            try:
                tax_id = str(seq_record.features[feat].qualifiers['db_xref'][0])
                organism = str(seq_record.features[feat].qualifiers['organism'][0])
                taxonomy_array.at[accession,'TaxonID'] = tax_id
                taxonomy_array.at[accession,'Organism'] = organism
            except KeyError:
                pass

taxonomy_array.TaxonID = taxonomy_array.TaxonID.str.replace(r'taxon:', '')

# Retrieve taxonomy information from Entrez

taxonomy_groups_to_test_lower = list(map(str.lower,taxonomy_groups_to_test))

current_accession_list = taxonomy_array.index
count = 0
length_of_list = len(current_accession_list)
for accession in current_accession_list:
    count += 1
    progress = str(count) + '/' + str(length_of_list)
    taxid = int(taxonomy_array.loc[accession,'TaxonID'])
    search = Entrez.efetch(id = taxid, db = "taxonomy", retmode = "xml")
    data = Entrez.read(search)
    lineage = {d['Rank']:d['ScientificName'] for d in data[0]['LineageEx'] if d['Rank'] in taxonomy_groups_to_test_lower}
    for i in range(len(taxonomy_groups_to_test)):
        if taxonomy_groups_to_test_lower[i] in lineage: taxonomy_array.at[accession,taxonomy_groups_to_test[i]] = lineage[taxonomy_groups_to_test_lower[i]]


edge_names = recip_TSS_plus_one_proteins_SUMS_BL['Edge_name']
internal_edge_names = [x for x in edge_names if "_" not in x]
query_file = internal_edge_names
tree_file = newick_tree_to_use

# Define functions for bringing consensus taxonomy information to internal edge

def lookup_by_names(tree):
    names = {}
    for clade in tree.find_clades():
        if clade.name:
            if clade.name in names:
                raise ValueError("Duplicate key: %s" % clade.name)
            names[clade.name] = clade
    return names

def children(tree):

    children_dict = {}
    for parent in tree.find_clades(terminal=False, order='level'):
        for child in parent.clades:
            if parent.name not in children_dict.keys():
                children_dict[parent.name] = [child.name]
            else:
                children_dict[parent.name].append(child.name)

    return children_dict

def get_leaves(tree):  # To get a list of all species in a given tree

    species = []
    for clade in tree.find_clades(terminal=True, order='level'):
        species.append(clade.name)

    return species

def process_query(query, leaves_list):

    decision_to_continue = True
    query = query.strip('\r')

    # Define the starting point based on the given input (Node vs. Edge)
    if '*' in query:  # Query is an Edge
        node1 = query.split('*')[0]
        node2 = query.split('*')[1]
        if clades_by_name[node1] in clades_by_name[node2]:
            starting_point = node1  # otherwise, descendent node is selected.
            remark = 'descendent'
            if starting_point in leaves_list:
                print ("The descendent node in query %s is a terminal leaf; no daughters exist." % query)
                decision_to_continue = False

        elif clades_by_name[node2] in clades_by_name[node1]:
            starting_point = node2  # otherwise, descendent node is selected.
            remark = 'descendent'
            if starting_point in leaves_list:
                print ("The descendent node in query %s is a terminal leaf; no daughters exist." % query)
                decision_to_continue = False
        else:
            print ("Are the input nodes connected? Please check your input.")
            sys.exit()
    else:  # Query is a Node
        starting_point = query
        remark = ''
        if starting_point in leaves_list:
            print ("The query node %s is a terminal leaf; no daughters exist." % query)
            decision_to_continue = False

    if decision_to_continue:

        if remark != '':
            outputname = 'Daughters_of_%s_in_Query_[%s]_in_%s_mode.txt' % (starting_point, query, remark)
        else:
            outputname = 'Daughters_of_%s_in_Query_[%s].txt' % (starting_point, query)

        paths_to_root = []

        for stop in my_tree.get_path(starting_point):
            paths_to_root.append(stop.name)  # Collect all nodes from the starting point up to the Root

        down_paths = []

        for leaf in leaves:
            leaf_path = []
            for stop in my_tree.get_path(leaf):
                leaf_path.append(stop.name)  # Collect all nodes from a given leaf up to the Root

            if starting_point in leaf_path:  # If the starting point is among the nodes collected above,
                for item in leaf_path:
                    if item not in paths_to_root:  # Record nodes that are between the leaf and the starting node
                        down_paths.append(item)

                        
        list_of_items = []
        for item in sorted(list(set(down_paths))):
            if item in children_dict.keys():
                pass
            else:
                list_of_items.append(item)
                number_of_daughters_X = len(list_of_items)
                
        #print ("Processing query: %s is completed." % query)
        return list_of_items,number_of_daughters_X

# Prepare taxonomy-related dataframes

taxonomy_daughters_of_internal_edge = pd.DataFrame(index = internal_edge_names, columns = taxonomy_columns_B)
taxonomy_daughters_of_internal_edge['Number_of_daughters_of_descendent_node'] = 0

# Read the input tree

my_tree = Phylo.read(tree_file, 'newick')

# Get all terminal leaves given in the input tree

leaves = get_leaves(my_tree)
clades_by_name = lookup_by_names(my_tree)

# Get the children of all nodes in the input tree

children_dict = children(my_tree)

for query in internal_edge_names:

    list_of_terminal_nodes, number_of_daughters = process_query(query, leaves)
    taxonomy_daughters_of_internal_edge.at[query,'Number_of_daughters_of_descendent_node'] = number_of_daughters
    working_taxonomy_dataframe = pd.DataFrame(list_of_terminal_nodes,columns=['Sequence_name'])
    working_taxonomy_dataframe = working_taxonomy_dataframe.set_index('Sequence_name').join(taxonomy_array.set_index('Sequence_name'))
    
    for taxonomic_level in taxonomy_columns_B:
        counts_of_entries = working_taxonomy_dataframe[taxonomic_level].value_counts()
        #print(counts_of_entries)
        if len(counts_of_entries) == 1:
            decision = counts_of_entries.index[0]
        elif len(counts_of_entries) == 0:
            decision = ''
        else:
            decision = 'MIXED'
        taxonomy_daughters_of_internal_edge.at[query,taxonomic_level] = decision
        print(query,decision)


# Further processing of taxonomy dataframes in preparation for merge

taxonomy_daughters_of_internal_edge.index.rename('Edge_name', inplace=True)
taxonomy_daughters_of_internal_edge.reset_index(inplace=True)
taxonomy_daughters_of_internal_edge[['Ancestral_node','Descendant_node_or_species']] = taxonomy_daughters_of_internal_edge.Edge_name.str.split("*",expand=True)
del taxonomy_daughters_of_internal_edge['Ancestral_node']

taxonomy_array['Number_of_daughters_of_descendent_node'] = 0
taxonomy_array.reset_index(drop=True,inplace=True)
taxonomy_array.drop('Organism',axis=1,inplace=True)
taxonomy_array.drop('TaxonID',axis=1,inplace=True)
taxonomy_array.rename({'Sequence_name':'Edge_name'},axis=1,inplace=True)
taxonomy_array['Descendant_node_or_species'] = taxonomy_array['Edge_name']

merge_on_total_taxonomy_array = ['Descendant_node_or_species','Edge_name','Number_of_daughters_of_descendent_node'] + taxonomy_groups_to_test
total_taxonomy_array = taxonomy_daughters_of_internal_edge.merge(taxonomy_array,on=merge_on_total_taxonomy_array,how='outer')

# Report everything about changes at selected positions

# To load the FASTA alignment file

def read_fasta(alignment):  # To read the FASTA alignment file
    aa_dict = {}
    with open(alignment, mode='r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            identifier = record.id
            sequence = record.seq
            aa_dict[identifier] = sequence
    return aa_dict

# To retrieve clades by name

def lookup_by_names(tree):
    names = {}
    for clade in tree.find_clades():
        if clade.name:
            if clade.name in names:
                raise ValueError("Duplicate key: %s" % clade.name)
            names[clade.name] = clade
    return names

# To retrieve all edges in the tree

def all_edges(tree):

    alledges = []
    for parent in tree.find_clades(terminal=False, order='level'):
        for child in parent.clades:
            alledges.append(parent.name + '*' + child.name)

    return alledges

def ISS_fluctuation_table_gap_threshold(PAGAN_tree,PAGAN_alignment,selected_protein,taxonomy_into_function):
    
    # Gap threshold

    record_x_toward_seq_dataframe = []
    sequence_records = []
    alignment_record_name_list = []

    for record in SeqIO.parse(PAGAN_alignment,"fasta"):
        alignment_record_name_list.append(record.name)
        record_x_toward_seq_dataframe = list(record.seq)
        record_x_toward_seq_dataframe_UPPER = [x.upper() for x in record_x_toward_seq_dataframe] 
        sequence_records.append(record_x_toward_seq_dataframe_UPPER)

    length_of_alignment = len(list(record.seq))
    depth_of_alignment = (len(alignment_record_name_list))

    accession_name_dataframe = pd.DataFrame(alignment_record_name_list, columns=['Accession'])
    sequence_dataframe = pd.DataFrame(sequence_records)

    gap_record = []
    sequence_columns = len(sequence_dataframe.axes[1])
    for i in range(sequence_columns):
        column_fractions_S = sequence_dataframe[i].value_counts(normalize='True')
        gap_fraction_column = 0
        for character, val in  column_fractions_S.iteritems():
            if character == '-':
                gap_fraction_column = val
        gap_record.append(gap_fraction_column)
    gap_fraction_DF = pd.DataFrame(np.nan, index=range(len(sequence_dataframe.axes[1])),columns=['Gap_fraction', 'Alignment_position'])
    gap_fraction_DF['Gap_fraction'] = gap_record
    gap_fraction_DF['Alignment_position'] = gap_fraction_DF.index.values + 1
    gap_fraction_DF_COPY = gap_fraction_DF.copy(deep=True)

    # Obtain alignment position series for later use

    all_alignment_positions = gap_fraction_DF_COPY.pop('Alignment_position')

    # Read the NWK PAGAN tree and record all edges

    my_tree = Phylo.read(PAGAN_tree, 'newick')
    edges = all_edges(my_tree)
    clades_dict = lookup_by_names(my_tree)

    # Read the PAGAN FASTA file and record character values for each position

    all_values = read_fasta(PAGAN_alignment)
    length_of_report_array = length_of_alignment * len(edges)
    dummy_array = np.empty((length_of_report_array,9),dtype=str)

    report_everything_about_selected_positions_results = pd.DataFrame(dummy_array,columns = ['Alignment_position','Edge_name','Edge_remark','Ancestral_node','Descendant_node_or_species','Ancestral_character','Descendant_character','Character_remark','Branch_length'])

    length_of_report_array_count = 0

    for i in range(length_of_alignment):
        query_label = str(i+1)
        query = i
        for edge in edges:
            parent_edge = edge.split('*')[0]
            child_edge = edge.split('*')[1]
            edge_len = str(clades_dict[child_edge].branch_length)
            parent_value = all_values[parent_edge][query]
            child_value = all_values[child_edge][query]
            if edge.count('_') == 0:
                edge_remark = 'Internal'
            elif edge.count('_') > 0:
                edge_remark = 'All'
            #else:
            #    print('Weird edge name', edge)
            #    sys.exit()
            if parent_value == child_value:
                value_remark = 'Conserved'
            else:
                value_remark = 'Fluctuating'

            to_append_report = ([query_label, edge, edge_remark, parent_edge, child_edge, parent_value, child_value,
                             value_remark, edge_len])
            report_everything_about_selected_positions_results.at[length_of_report_array_count]=to_append_report
            length_of_report_array_count += 1

    report_everything_about_selected_positions_results = report_everything_about_selected_positions_results[report_everything_about_selected_positions_results['Character_remark'] != 'Conserved']

    fluctuation_table = report_everything_about_selected_positions_results.copy()

    # Place gap metrics in fluctuation table with a join on 'Alignment position'

    fluctuation_table['Alignment_position'] = fluctuation_table['Alignment_position'].astype(int)
    fluctuation_table = fluctuation_table.set_index('Alignment_position').join(gap_fraction_DF.set_index('Alignment_position'))
    fluctuation_table = fluctuation_table.reset_index()

    # Remove entries based on gap fraction

    fluctuation_table_bool_series_remove = fluctuation_table['Gap_fraction'] > gap_cutoff
    fluctuation_table_bool_index_remove = fluctuation_table_bool_series_remove[fluctuation_table_bool_series_remove].index
    fluctuation_table = fluctuation_table.drop(fluctuation_table_bool_index_remove)

    # Remove entries based on ancestral or descendant character as non-standard protein

    fluctuation_table['Ancestral_character'] = fluctuation_table['Ancestral_character'].astype(str)
    fluctuation_table['Descendant_character'] = fluctuation_table['Descendant_character'].astype(str)
    fluctuation_table = fluctuation_table.loc[~fluctuation_table['Ancestral_character'].isin(['X','-'])]
    fluctuation_table = fluctuation_table.loc[~fluctuation_table['Descendant_character'].isin(['X','-'])]

    # Remove gap percent column - no longer needed

    del fluctuation_table['Gap_fraction']
    fluctuation_table = fluctuation_table.reset_index(drop=True)
    fluctuation_table.head()

    # Count fluctuations at each position across internal and external edges

    counts_relevant_TSS = fluctuation_table['Alignment_position'].value_counts()
    counts_relevant_TSS = counts_relevant_TSS.to_frame()

    counts_relevant_TSS.reset_index(level=0, inplace=True)
    counts_relevant_TSS.columns = ['Alignment_position','Total_substitution_score']

    all_alignment_positions_DF = all_alignment_positions.to_frame()

    # Prepare a fluctuation table across internal and external edges (TSS) (below gap threshold and no 'X' or '-')

    TSS_fluctuation_table = fluctuation_table.merge(counts_relevant_TSS,on='Alignment_position',how='left')
    
    # Total substitution score table by alignment positions
    
    total_substitution_score = all_alignment_positions_DF.merge(counts_relevant_TSS,on='Alignment_position',how='left')
    total_substitution_score['Total_substitution_score'] = total_substitution_score['Total_substitution_score'].fillna(0)
    total_substitution_score['Total_substitution_score'].iloc[gap_fraction_DF['Gap_fraction'] > gap_cutoff] = 'HIGH_GAP: >'+str(gap_cutoff*100)+'%'
    #total_substitution_score = total_substitution_score.reset_index()
   
    total_substitution_score_gap_removed = total_substitution_score[total_substitution_score['Total_substitution_score'] != 'HIGH_GAP: >'+str(gap_cutoff*100)+'%']
    
    # Remove external tree edges

    external_edges_bool = fluctuation_table['Edge_name'].str.contains('_') # '-' is the character found in external edges
    external_edges_index_remove = external_edges_bool[external_edges_bool].index
    fluctuation_table = fluctuation_table.drop(external_edges_index_remove)
    fluctuation_table = fluctuation_table.reset_index(drop=True)

    # Count fluctuations at each position to generate internal substitution score

    counts_relevant_ISS = fluctuation_table['Alignment_position'].value_counts()
    counts_relevant_ISS = counts_relevant_ISS.to_frame()

    counts_relevant_ISS.reset_index(level=0, inplace=True)
    counts_relevant_ISS.columns = ['Alignment_position','Internal_substitution_score']
    
    # Prepare fluctuation table for substitutions along internal edges that contains ISS values

    ISS_fluctuation_table = fluctuation_table.merge(counts_relevant_ISS,on='Alignment_position',how='left')
    
    # Generate an ISS and alignment table, first saved with gap notations, then remove positions above gap threshold
    
    internal_substitution_score = all_alignment_positions_DF.merge(counts_relevant_ISS,on='Alignment_position',how='left')
    internal_substitution_score['Internal_substitution_score'] = internal_substitution_score['Internal_substitution_score'].fillna(0)
    internal_substitution_score['Internal_substitution_score'].iloc[gap_fraction_DF['Gap_fraction'] > gap_cutoff] = 'HIGH_GAP: >'+str(gap_cutoff*100)+'%'
    
    internal_substitution_score_gap_removed = internal_substitution_score[internal_substitution_score['Internal_substitution_score'] != 'HIGH_GAP: >'+str(gap_cutoff*100)+'%']
    
    # Generate tables for conservation versus alignment position
    
    conservation_table = internal_substitution_score.merge(total_substitution_score,on='Alignment_position',how='left')
    #conservation_table.reset_index(drop=True, inplace=True)
    conservation_table.to_csv(file_prefix + 'conservation_table_TSS_and_ISS_'+selected_protein+'.csv', index=False)
    
    # Add taxonomy data to fluctuation_tables
    
    #ISS_fluctuation_table.to_csv('ISS_fluctuation_table_'+selected_protein+'.tabular',sep ='\t', index=False)
    
    taxonomy_into_function_RC = taxonomy_into_function.copy()
    del taxonomy_into_function_RC['Edge_name']

    TSS_fluctuation_table_with_taxonomy = TSS_fluctuation_table.merge(taxonomy_into_function_RC,on='Descendant_node_or_species',how='left')
    ISS_fluctuation_table_with_taxonomy = ISS_fluctuation_table.merge(taxonomy_into_function_RC,on='Descendant_node_or_species',how='left')
    
    return TSS_fluctuation_table_with_taxonomy,ISS_fluctuation_table_with_taxonomy,total_substitution_score_gap_removed


# Save total taxonomy array

total_taxonomy_array.to_csv(file_prefix + 'total_taxonomy_array.csv',index=False)

# Cycle through selected proteins

ALL_PROTEINS_returned_TSS_fluctuation = pd.DataFrame(columns=['Alignment_position','Protein','Edge_name','Edge_remark','Ancestral_node',\
                                                              'Descendant_node_or_species','Ancestral_character','Descendant_character','Character_remark',\
                                                              'Branch_length','Total_substitution_score'] + taxonomy_groups_to_test + ['Number_of_daughters_of_descendent_node','recip_TSS_plus_one'])
ALL_PROTEINS_returned_ISS_fluctuation = pd.DataFrame(columns=['Alignment_position','Protein','Edge_name','Edge_remark','Ancestral_node',\
                                                              'Descendant_node_or_species','Ancestral_character','Descendant_character','Character_remark',\
                                                              'Branch_length','Internal_substitution_score'] + taxonomy_groups_to_test + ['Number_of_daughters_of_descendent_node'])

list_of_all_ISS_values = []
list_of_all_TSS_values = []
recip_ISS_plus_one_columns = []
recip_TSS_plus_one_columns = []
print('Proteins analyzed: ')
for protein in chosen_genes:
    PAGAN_sequence_to_test = file_prefix + protein + '_prot_PAGAN_'+ selected_accession + '.fasta'
    PAGAN_tree_to_test = file_prefix + protein + '_prot_PAGAN.anctree'
    fluct_frame = 'fluctframe_' + protein
    
    TSS_position_frame = 'TSSpos_' + protein
    returned_TSS_fluctuation, returned_ISS_fluctuation, returned_TSS_column = ISS_fluctuation_table_gap_threshold(PAGAN_tree_to_test,PAGAN_sequence_to_test,protein,total_taxonomy_array) ##
    returned_TSS_fluctuation['recip_TSS_plus_one'] = 1/(returned_TSS_fluctuation['Total_substitution_score']+1)
    sum_recip_TSS_plus_one = returned_TSS_fluctuation.groupby(['Edge_name']).recip_TSS_plus_one.sum().reset_index()
    sum_recip_TSS_plus_one.rename(columns = {'recip_TSS_plus_one':protein+'_recip_TSS_plus_one'}, inplace = True)
    recip_TSS_plus_one_columns.append(protein+'_recip_TSS_plus_one')
    sum_for_merge = sum_recip_TSS_plus_one.copy()
    recip_TSS_plus_one_proteins_SUMS_BL = recip_TSS_plus_one_proteins_SUMS_BL.merge(sum_for_merge, on='Edge_name', how='outer')
    globals()[TSS_position_frame] = returned_TSS_column
    TSS_pos_list = globals()[TSS_position_frame]['Total_substitution_score'].values.tolist()
    list_of_all_TSS_values = list_of_all_TSS_values + TSS_pos_list
    
    returned_TSS_fluctuation['Protein'] = protein
    returned_ISS_fluctuation['Protein'] = protein
    ALL_PROTEINS_returned_TSS_fluctuation = ALL_PROTEINS_returned_TSS_fluctuation.append(returned_TSS_fluctuation, ignore_index=False)
    ALL_PROTEINS_returned_ISS_fluctuation = ALL_PROTEINS_returned_ISS_fluctuation.append(returned_ISS_fluctuation, ignore_index=False)
    print(protein+': '+str(len(returned_TSS_fluctuation)), 'substitutions along internal and external edges at analyzed positions.')
    

# Save total and internal fluctuation tables that contain all analyzed proteins

ALL_PROTEINS_returned_TSS_fluctuation.drop(['Character_remark','recip_TSS_plus_one'],axis=1,inplace=True)
ALL_PROTEINS_returned_ISS_fluctuation.drop(['Character_remark'],axis=1,inplace=True)

ALL_PROTEINS_returned_TSS_fluctuation = ALL_PROTEINS_returned_TSS_fluctuation.sort_values(by='Total_substitution_score', ascending=True)
ALL_PROTEINS_returned_ISS_fluctuation = ALL_PROTEINS_returned_ISS_fluctuation.sort_values(by='Internal_substitution_score', ascending=True)

ALL_PROTEINS_returned_TSS_fluctuation.to_csv(file_prefix + 'TSS_fluctuation_table_all_proteins.csv', index=False)
ALL_PROTEINS_returned_ISS_fluctuation.to_csv(file_prefix + 'ISS_fluctuation_table_all_proteins.csv', index=False)

# Graph ISS values across all tested positions and save TSS distributions

sorted_list_of_all_TSS_values = sorted(list_of_all_TSS_values)
sorted_TSS_values_S = pd.Series(sorted_list_of_all_TSS_values,name='Instances_of_TSS')
# sorted_TSS_values_S.plot(xlabel = 'Protein site analyzed (ordered with increasing TSS value)',ylabel = 'TSS')

sorted_TSS_values_S.to_csv(file_prefix + 'all_site_TSS_across_analysis.csv',index=False)

value_counts_pos_TSS = sorted_TSS_values_S.value_counts(sort=True, ascending=True)
value_counts_pos_TSS.index.name = 'TSS'

value_counts_pos_TSS.to_csv(file_prefix + 'frequency_of_TSS.csv')

del recip_TSS_plus_one_proteins_SUMS_BL['Distance_from_root']
recip_TSS_plus_one_proteins_SUMS_BL.fillna(0, inplace=True)

# All protein scores sum

recip_TSS_plus_one_proteins_SUMS_BL['Sum-of-protein-scores_recip_TSS_plus_one']=recip_TSS_plus_one_proteins_SUMS_BL.loc[:,recip_TSS_plus_one_columns].sum(axis=1)

# Complex-specific sums for OXPHOS complexes (unquote to allow)


recip_TSS_plus_one_proteins_SUMS_BL['Complex-I_recip_TSS_plus_one'] = recip_TSS_plus_one_proteins_SUMS_BL['ND1_recip_TSS_plus_one'] + \
    recip_TSS_plus_one_proteins_SUMS_BL['ND2_recip_TSS_plus_one'] + \
    recip_TSS_plus_one_proteins_SUMS_BL['ND3_recip_TSS_plus_one'] + \
    recip_TSS_plus_one_proteins_SUMS_BL['ND4_recip_TSS_plus_one'] + \
    recip_TSS_plus_one_proteins_SUMS_BL['ND4L_recip_TSS_plus_one'] + \
    recip_TSS_plus_one_proteins_SUMS_BL['ND5_recip_TSS_plus_one'] + \
    recip_TSS_plus_one_proteins_SUMS_BL['ND6_recip_TSS_plus_one']

recip_TSS_plus_one_proteins_SUMS_BL['Complex-III_recip_TSS_plus_one'] = recip_TSS_plus_one_proteins_SUMS_BL['CYTB_recip_TSS_plus_one']

recip_TSS_plus_one_proteins_SUMS_BL['Complex-IV_recip_TSS_plus_one'] = recip_TSS_plus_one_proteins_SUMS_BL['COX1_recip_TSS_plus_one'] + \
    recip_TSS_plus_one_proteins_SUMS_BL['COX2_recip_TSS_plus_one'] + \
    recip_TSS_plus_one_proteins_SUMS_BL['COX3_recip_TSS_plus_one']

recip_TSS_plus_one_proteins_SUMS_BL['Complex-V_recip_TSS_plus_one'] = recip_TSS_plus_one_proteins_SUMS_BL['ATP8_recip_TSS_plus_one'] + \
    recip_TSS_plus_one_proteins_SUMS_BL['ATP6_recip_TSS_plus_one']

recip_TSS_plus_one_proteins_SUMS_BL = recip_TSS_plus_one_proteins_SUMS_BL.sort_values(by='Sum-of-protein-scores_recip_TSS_plus_one', ascending=False)
recip_TSS_plus_one_proteins_SUMS_BL.to_csv(file_prefix + 'tree_edge_and_protein_sums_recip_TSS_plus_one.csv',index=False)

# Remove those entries with maximum branch lengths (typically 0.2 when using PAGAN) - can change as desired - prevents linear regression using samples with highly aberrant branch lengths

maximum_branch_length = recip_TSS_plus_one_proteins_SUMS_BL['Branch_length'].max()
recip_TSS_plus_one_proteins_SUMS_BL = recip_TSS_plus_one_proteins_SUMS_BL[recip_TSS_plus_one_proteins_SUMS_BL['Branch_length'] != maximum_branch_length]

# Scaling 0 to 1


columns_for_scaling = recip_TSS_plus_one_proteins_SUMS_BL.columns[recip_TSS_plus_one_proteins_SUMS_BL.columns.str.contains(pat = '_recip_TSS_plus_one')].tolist()
scaler = MinMaxScaler()
recip_TSS_plus_one_proteins_SUMS_BL['Branch_length_scaled'] = scaler.fit_transform(recip_TSS_plus_one_proteins_SUMS_BL['Branch_length'].values.reshape(-1,1))
sep = '_'

for column in columns_for_scaling:
    stem = column.split(sep, 1)[0]
    recip_TSS_plus_one_proteins_SUMS_BL[stem+'_scaled'] = scaler.fit_transform(recip_TSS_plus_one_proteins_SUMS_BL[column].values.reshape(-1,1))
    
recip_TSS_plus_one_proteins_SUMS_BL['Sum-of-protein-scores_scaled'] = scaler.fit_transform(recip_TSS_plus_one_proteins_SUMS_BL['Sum-of-protein-scores_recip_TSS_plus_one'].values.reshape(-1,1))

# Linear regression

columns_for_regression = columns_for_scaling

for column in columns_for_regression:
    
    X = recip_TSS_plus_one_proteins_SUMS_BL['Branch_length_scaled'].to_numpy()
    stem = column.split(sep, 1)[0]
    y = recip_TSS_plus_one_proteins_SUMS_BL[stem + '_scaled'].to_numpy()
    
    idx = np.isfinite(X) & np.isfinite(y)
    coeffs = np.polyfit(X[idx], y[idx], 1)
    
    coeff_x_power_1 = coeffs[0]
    coeff_0 = coeffs[1]
    
    recip_TSS_plus_one_proteins_SUMS_BL[stem + '_expected'] = coeff_x_power_1 * recip_TSS_plus_one_proteins_SUMS_BL['Branch_length_scaled'] + coeff_0


# Residuals

for column in columns_for_regression:
    stem = column.split(sep, 1)[0]
    recip_TSS_plus_one_proteins_SUMS_BL[stem + '_residual'] = recip_TSS_plus_one_proteins_SUMS_BL[stem + '_scaled'] - recip_TSS_plus_one_proteins_SUMS_BL[stem + '_expected']

# Arctangent

for column in columns_for_regression:
    stem = column.split(sep, 1)[0]
    recip_TSS_plus_one_proteins_SUMS_BL[stem + '_arctan'] = np.arctan(recip_TSS_plus_one_proteins_SUMS_BL[stem + '_residual']/recip_TSS_plus_one_proteins_SUMS_BL['Branch_length_scaled'])

# Rank percentile:

for column in columns_for_regression:
    stem = column.split(sep, 1)[0]
    recip_TSS_plus_one_proteins_SUMS_BL[stem + '_rankpercent'] = recip_TSS_plus_one_proteins_SUMS_BL[stem + '_arctan'].rank(pct=True, ascending=False)

# Calculate arctan values for tree edges

del total_taxonomy_array['Edge_name']

edge_and_arctan_columns = recip_TSS_plus_one_proteins_SUMS_BL.columns[recip_TSS_plus_one_proteins_SUMS_BL.columns.str.contains(pat = '_arctan')].tolist()
columns_selected_for_arctan_DF = ['Edge_name'] + edge_and_arctan_columns
edge_and_arctan = recip_TSS_plus_one_proteins_SUMS_BL[columns_selected_for_arctan_DF]
edge_and_arctan[['Ancestral_node','Descendant_node_or_species']] = edge_and_arctan.Edge_name.str.split("*",expand=True)
edge_and_arctan = edge_and_arctan.merge(total_taxonomy_array,on=['Descendant_node_or_species'],how='left')
edge_and_arctan = edge_and_arctan.sort_values(by='Sum-of-protein-scores_arctan', ascending=False)
edge_and_arctan.to_csv(file_prefix + 'tree_edge_and_protein_arctangents.csv',index=False)

# Calculate rankpercent values for tree edges

edge_and_rankperc_columns = recip_TSS_plus_one_proteins_SUMS_BL.columns[recip_TSS_plus_one_proteins_SUMS_BL.columns.str.contains(pat = '_rankpercent')].tolist()
columns_selected_for_rankperc_DF = ['Edge_name'] + edge_and_rankperc_columns
edge_and_rankperc = recip_TSS_plus_one_proteins_SUMS_BL[columns_selected_for_rankperc_DF]
edge_and_rankperc[['Ancestral_node','Descendant_node_or_species']] = edge_and_rankperc.Edge_name.str.split("*",expand=True)
edge_and_rankperc = edge_and_rankperc.merge(total_taxonomy_array,on=['Descendant_node_or_species'],how='left')
edge_and_rankperc = edge_and_rankperc.sort_values(by='Sum-of-protein-scores_rankpercent', ascending=True)
edge_and_rankperc.to_csv(file_prefix + 'tree_edge_and_protein_rankpercentiles.csv',index=False)

# Confidence levels (with assistance from https://www.geeksforgeeks.org/how-to-plot-a-confidence-interval-in-python/)

list_of_results = []

for group in taxonomy_groups_to_test:
    clade_S = pd.Series(edge_and_arctan[group])
    clade_S = clade_S[~clade_S.isnull()]
    clades = set(clade_S.to_list())
    for clade in clades:
        for column in arctan_columns_to_assess:
            copy_of_arctangent_DF = edge_and_arctan.copy()
            selected_rows_by_spec_group = copy_of_arctangent_DF[copy_of_arctangent_DF[group] == clade] # select rows by that particular clade
            sum_prot_scores_arc_NP = selected_rows_by_spec_group[column].to_numpy()
            median_of_clade = np.median(sum_prot_scores_arc_NP)
            n_iterations = 1000 # here k=no. of bootstrapped samples 
            n_size = int(len(sum_prot_scores_arc_NP)) 
            if clade == 'nan': pass
            # run bootstrap 
            medians = list() 
            for i in range(n_iterations): 
               s = resample(sum_prot_scores_arc_NP, n_samples=n_size); 
               m = np.median(s); 
               medians.append(m) 

            p = ((1.0-alpha)/2.0) * 100
            lower =  np.percentile(medians, p) 
            p = (alpha+((1.0-alpha)/2.0)) * 100
            upper =  np.percentile(medians, p) 
            range_CI = upper - lower

            data_to_add = group,clade,column,median_of_clade,lower,upper,n_size,range_CI
            list_of_results.append(data_to_add)
alpha_to_percent = str(int(alpha * 100))       
dataframe_of_results = pd.DataFrame(list_of_results,columns = ['Taxonomic_Level','Clade','Arctan_calculation','Median_of_Clade',alpha_to_percent+'%_CI_Lower',alpha_to_percent+'%_CI_Upper','Number_of_Samples','CI_Range'])

for column in arctan_columns_to_assess:
    copy_of_dataframe_of_results = dataframe_of_results.copy()
    copy_of_dataframe_of_results = copy_of_dataframe_of_results[copy_of_dataframe_of_results['Arctan_calculation'] == column]
    copy_of_dataframe_of_results.to_csv(file_prefix + column + '.csv',index=False)

# Print lists of values by taxonomic level and clade

for group in taxonomy_groups_to_test:
    clade_S = pd.Series(edge_and_arctan[group])
    clade_S = clade_S[~clade_S.isnull()]
    clades = set(clade_S.to_list())
    for column in arctan_columns_to_assess:
        column_names = []
        d = []
        for clade in clades:
            copy_of_arctangent_DF = edge_and_arctan.copy()
            selected_rows_by_spec_group = copy_of_arctangent_DF[copy_of_arctangent_DF[group] == clade] # select rows by that particular clade
            sum_prot_scores_arc_list = selected_rows_by_spec_group[column].to_list()
            globals()[clade] = sum_prot_scores_arc_list
            column_names.append(clade)
            d.append(globals()[clade])
        export_data = zip_longest(*d, fillvalue = '')
        with open(file_prefix + group + '_' + column + '_for_prism.csv', 'w', encoding="ISO-8859-1", newline='') as myfile:
            wr = csv.writer(myfile)
            wr.writerow((column_names))
            wr.writerows(export_data)
        myfile.close()

