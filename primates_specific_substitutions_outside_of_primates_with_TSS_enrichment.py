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

# Read full substitution table and generate a column fusing gene and substitution

total_fluctuations_TSS = pd.read_csv('mammal_and_Ap_mtDNA_PYTEST_FEB_21_TSS_fluctuation_table_all_proteins.csv')
total_fluctuations_TSS['Substitution'] = total_fluctuations_TSS['Protein'].astype(str) + '_' + total_fluctuations_TSS['Ancestral_character'].astype(str) + total_fluctuations_TSS['Alignment_position'].astype(str) + total_fluctuations_TSS['Descendant_character'].astype(str)

# Generate a dataframe containing those substitutions found in primates

primates_substitutions = total_fluctuations_TSS[(total_fluctuations_TSS['Order'] == 'Primates')]

# Generate a dataframe containing those substitutions not found in primates

not_primates_substitutions = total_fluctuations_TSS[(total_fluctuations_TSS['Order'] != 'Primates')]

# Generate a list of substitutions from the primates substitution dataset

primates_substitutions_set = list(set(primates_substitutions['Substitution'].to_list()))

# Select those substitutions outside of primates that were matched (ancestral character > descendent character) inside the primates dataset

not_primates_subsitutions_primates_limited = not_primates_substitutions[not_primates_substitutions['Substitution'].isin(primates_substitutions_set)]
not_primates_subsitutions_primates_limited.to_csv('not_primates_subsitutions_matching_those_in_primates_limited.csv',index=False)

# Filter this dataset based upon TSS limit

not_primates_subsitutions_primates_limited_100 = not_primates_subsitutions_primates_limited[not_primates_subsitutions_primates_limited['Total_substitution_score'] <= 100]
not_primates_subsitutions_primates_limited_50 = not_primates_subsitutions_primates_limited[not_primates_subsitutions_primates_limited['Total_substitution_score'] <= 50]
not_primates_subsitutions_primates_limited_20 = not_primates_subsitutions_primates_limited[not_primates_subsitutions_primates_limited['Total_substitution_score'] <= 20]
not_primates_subsitutions_primates_limited_10 = not_primates_subsitutions_primates_limited[not_primates_subsitutions_primates_limited['Total_substitution_score'] <= 10]
not_primates_subsitutions_primates_limited_5 = not_primates_subsitutions_primates_limited[not_primates_subsitutions_primates_limited['Total_substitution_score'] <= 5]

# Determine fraction of these substitutions found in each order at different TSS limits

not_primates_subsitutions_primates_limited_ALL_ordercounts = not_primates_subsitutions_primates_limited['Order'].value_counts(normalize=True)
not_primates_subsitutions_primates_limited_100_ordercounts = not_primates_subsitutions_primates_limited_100['Order'].value_counts(normalize=True)
not_primates_subsitutions_primates_limited_50_ordercounts = not_primates_subsitutions_primates_limited_50['Order'].value_counts(normalize=True)
not_primates_subsitutions_primates_limited_20_ordercounts = not_primates_subsitutions_primates_limited_20['Order'].value_counts(normalize=True)
not_primates_subsitutions_primates_limited_10_ordercounts = not_primates_subsitutions_primates_limited_10['Order'].value_counts(normalize=True)
not_primates_subsitutions_primates_limited_5_ordercounts = not_primates_subsitutions_primates_limited_5['Order'].value_counts(normalize=True)

# Determine fraction of these substitutions found in each family at different TSS limits

not_primates_subsitutions_primates_limited_ALL_familycounts = not_primates_subsitutions_primates_limited['Family'].value_counts(normalize=True)
not_primates_subsitutions_primates_limited_100_familycounts = not_primates_subsitutions_primates_limited_100['Family'].value_counts(normalize=True)
not_primates_subsitutions_primates_limited_50_familycounts = not_primates_subsitutions_primates_limited_50['Family'].value_counts(normalize=True)
not_primates_subsitutions_primates_limited_20_familycounts = not_primates_subsitutions_primates_limited_20['Family'].value_counts(normalize=True)
not_primates_subsitutions_primates_limited_10_familycounts = not_primates_subsitutions_primates_limited_10['Family'].value_counts(normalize=True)
not_primates_subsitutions_primates_limited_5_familycounts = not_primates_subsitutions_primates_limited_5['Family'].value_counts(normalize=True)

# Determine fraction of these substitutions found in each genus at different TSS limits

not_primates_subsitutions_primates_limited_ALL_genuscounts = not_primates_subsitutions_primates_limited['Genus'].value_counts(normalize=True)
not_primates_subsitutions_primates_limited_100_genuscounts = not_primates_subsitutions_primates_limited_100['Genus'].value_counts(normalize=True)
not_primates_subsitutions_primates_limited_50_genuscounts = not_primates_subsitutions_primates_limited_50['Genus'].value_counts(normalize=True)
not_primates_subsitutions_primates_limited_20_genuscounts = not_primates_subsitutions_primates_limited_20['Genus'].value_counts(normalize=True)
not_primates_subsitutions_primates_limited_10_genuscounts = not_primates_subsitutions_primates_limited_10['Genus'].value_counts(normalize=True)
not_primates_subsitutions_primates_limited_5_genuscounts = not_primates_subsitutions_primates_limited_5['Genus'].value_counts(normalize=True)

 # Build dataframes to hold fractional data for orders

fraction_across_order = pd.DataFrame({'All': not_primates_subsitutions_primates_limited_ALL_ordercounts,\
                                      '100': not_primates_subsitutions_primates_limited_100_ordercounts,\
                                      '50': not_primates_subsitutions_primates_limited_50_ordercounts,\
                                      '20':not_primates_subsitutions_primates_limited_20_ordercounts,\
                                      '10':not_primates_subsitutions_primates_limited_10_ordercounts,\
                                      '5':not_primates_subsitutions_primates_limited_5_ordercounts})
fraction_across_order = fraction_across_order.fillna(0)
fraction_across_order = fraction_across_order.drop('MIXED')
fraction_across_order = fraction_across_order.drop('Squamata')

 # Calculate log2 fold change for TSS-limited data versus data unlimited by TSS for orders

fraction_across_order['Log2_Fold_5_vs_All'] = np.log2(fraction_across_order['5']/fraction_across_order['All'])
fraction_across_order['Log2_Fold_10_vs_All'] = np.log2(fraction_across_order['10']/fraction_across_order['All'])
fraction_across_order['Log2_Fold_20_vs_All'] = np.log2(fraction_across_order['20']/fraction_across_order['All'])
fraction_across_order['Log2_Fold_50_vs_All'] = np.log2(fraction_across_order['50']/fraction_across_order['All'])
fraction_across_order['Log2_Fold_100_vs_All'] = np.log2(fraction_across_order['100']/fraction_across_order['All'])
fraction_across_order = fraction_across_order.sort_values('Log2_Fold_5_vs_All')

 # Build dataframes to hold fractional data for families

fraction_across_family = pd.DataFrame({'All': not_primates_subsitutions_primates_limited_ALL_familycounts,\
                                      '100': not_primates_subsitutions_primates_limited_100_familycounts,\
                                      '50': not_primates_subsitutions_primates_limited_50_familycounts,\
                                      '20':not_primates_subsitutions_primates_limited_20_familycounts,\
                                      '10':not_primates_subsitutions_primates_limited_10_familycounts,\
                                      '5':not_primates_subsitutions_primates_limited_5_familycounts})
fraction_across_family = fraction_across_family.fillna(0)
fraction_across_family = fraction_across_family.drop('MIXED')
fraction_across_family = fraction_across_family.drop('Dactyloidae')

 # Calculate log2 fold change for TSS-limited data versus data unlimited by TSS for families

fraction_across_family['Log2_Fold_5_vs_All'] = np.log2(fraction_across_family['5']/fraction_across_family['All'])
fraction_across_family['Log2_Fold_10_vs_All'] = np.log2(fraction_across_family['10']/fraction_across_family['All'])
fraction_across_family['Log2_Fold_20_vs_All'] = np.log2(fraction_across_family['20']/fraction_across_family['All'])
fraction_across_family['Log2_Fold_50_vs_All'] = np.log2(fraction_across_family['50']/fraction_across_family['All'])
fraction_across_family['Log2_Fold_100_vs_All'] = np.log2(fraction_across_family['100']/fraction_across_family['All'])
fraction_across_family = fraction_across_family.sort_values('Log2_Fold_5_vs_All')

 # Build dataframes to hold fractional data for genera

fraction_across_genus = pd.DataFrame({'All': not_primates_subsitutions_primates_limited_ALL_genuscounts,\
                                      '100': not_primates_subsitutions_primates_limited_100_genuscounts,\
                                      '50': not_primates_subsitutions_primates_limited_50_genuscounts,\
                                      '20':not_primates_subsitutions_primates_limited_20_genuscounts,\
                                      '10':not_primates_subsitutions_primates_limited_10_genuscounts,\
                                      '5':not_primates_subsitutions_primates_limited_5_genuscounts})
fraction_across_genus = fraction_across_genus.fillna(0)
fraction_across_genus = fraction_across_genus.drop('MIXED')
fraction_across_genus = fraction_across_genus.drop('Anolis')

# Calculate log2 fold change for TSS-limited data versus data unlimited by TSS for genera

fraction_across_genus['Log2_Fold_5_vs_All'] = np.log2(fraction_across_genus['5']/fraction_across_genus['All'])
fraction_across_genus['Log2_Fold_10_vs_All'] = np.log2(fraction_across_genus['10']/fraction_across_genus['All'])
fraction_across_genus['Log2_Fold_20_vs_All'] = np.log2(fraction_across_genus['20']/fraction_across_genus['All'])
fraction_across_genus['Log2_Fold_50_vs_All'] = np.log2(fraction_across_genus['50']/fraction_across_genus['All'])
fraction_across_genus['Log2_Fold_100_vs_All'] = np.log2(fraction_across_genus['100']/fraction_across_genus['All'])
fraction_across_genus = fraction_across_genus.sort_values('Log2_Fold_5_vs_All')

# Save to CSV

fraction_across_order.to_csv('primates_Substitution_Fraction_Across_Order_Log2_Fold_All_100_50_20_10_5.csv')
fraction_across_family.to_csv('primates_Substitution_Fraction_Across_Family_Log2_Fold_All_100_50_20_10_5.csv')
fraction_across_genus.to_csv('primates_Substitution_Fraction_Across_Genus_Log2_Fold_All_100_50_20_10_5.csv')

