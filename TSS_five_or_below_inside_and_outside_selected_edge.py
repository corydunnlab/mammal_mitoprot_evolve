#!/usr/bin/env python
# coding: utf-8 

# Funding received from the European Research Council, the Sigrid Jusélius Foundation, and the Academy of Finland contributed to the development of this software.
# Author: Cory Dunn
# Institution: University of Helsinki
# Author Email: cory.dunn@helsinki.fi
# License: GPLv3

# Import dependencies

import pandas as pd
import numpy as np

selected_edge = '#871#*#870#' # <----- Selected edge
node1,node2 = selected_edge.split('*')
TSS_limit = 5 # <----- Selected TSS limit

total_fluctuations_TSS = pd.read_csv('mammal_and_Ap_mtDNA_PYTEST_FEB_21_TSS_fluctuation_table_all_proteins.csv')
total_fluctuations_TSS['Substitution'] = total_fluctuations_TSS['Protein'].astype(str) + '_' + total_fluctuations_TSS['Ancestral_character'].astype(str) + total_fluctuations_TSS['Alignment_position'].astype(str) + total_fluctuations_TSS['Descendant_character'].astype(str)
total_fluctuations_TSS['Reverse_substitution'] = total_fluctuations_TSS['Protein'].astype(str) + '_' + total_fluctuations_TSS['Descendant_character'].astype(str) + total_fluctuations_TSS['Alignment_position'].astype(str) + total_fluctuations_TSS['Ancestral_character'].astype(str)

# Extract all substitutions within edge that are at positions at or below TSS limit and save to file

select_edge_substitutions = total_fluctuations_TSS.loc[(total_fluctuations_TSS['Edge_name'] == selected_edge) \
                                                          & (total_fluctuations_TSS['Total_substitution_score'] <= TSS_limit)]

select_edge_substitutions_set = list(set(select_edge_substitutions['Substitution'].to_list()))

matching_all_substitutions_within_this_edge = total_fluctuations_TSS.loc[(total_fluctuations_TSS['Substitution'].isin(select_edge_substitutions_set)) & (total_fluctuations_TSS['Edge_name'] == selected_edge)]

matching_all_substitutions_within_this_edge.to_csv('Substitutions_matching_'+selected_edge+'_TSS_lessthanorequalto_'+str(TSS_limit)+'_inside_edge.csv',index=False)

# Find all reversions across all edges matching substitutions below selected TSS in selected edge

select_reversion_set = list(set(select_edge_substitutions['Reverse_substitution'].to_list()))

matching_all_reversions = total_fluctuations_TSS.loc[(total_fluctuations_TSS['Substitution'].isin(select_reversion_set))]

matching_all_reversions.to_csv('Potential_reversions_of_'+selected_edge+'_TSS_lessthanorequalto_'+str(TSS_limit)+'_across_all_edges.csv',index=False)

# Extract all substitutions across all edges that match selected edge and TSS (but not within edge)

matching_substitutions_across_all_edges_but_this_one = total_fluctuations_TSS.loc[(total_fluctuations_TSS['Substitution'].isin(select_edge_substitutions_set)) & (total_fluctuations_TSS['Edge_name'] != selected_edge)]
matching_substitutions_across_all_edges_but_this_one.to_csv('Substitutions_matching_'+selected_edge+'_TSS_lessthanorequalto_'+str(TSS_limit)+'_outside_selected_edge.csv',index=False)

