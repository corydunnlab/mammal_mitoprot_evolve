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

# Load total substitution table

total_fluctuation_table = pd.read_csv('/Users/corydunn/Dropbox/Lab/Current_Lab_Members/CDUNN/University_of_Helsinki/CDD_Mammal_Selection_JAN_25_2021/python_testing_ground/mammal_and_Ap_mtDNA_PYTEST_FEB_21_TSS_fluctuation_table_all_proteins.csv')

# Break down total substitution table by TSS limit

total_fluctuation_table_five = total_fluctuation_table[total_fluctuation_table['Total_substitution_score'] <= 5]
total_fluctuation_table_ten = total_fluctuation_table[total_fluctuation_table['Total_substitution_score'] <= 10]
total_fluctuation_table_twenty = total_fluctuation_table[total_fluctuation_table['Total_substitution_score'] <= 20]
total_fluctuation_table_fifty = total_fluctuation_table[total_fluctuation_table['Total_substitution_score'] <= 50]
total_fluctuation_table_hundred = total_fluctuation_table[total_fluctuation_table['Total_substitution_score'] <= 100]

# Get the fraction of substitutions in each order at decreasing TSS

fraction_five_order = total_fluctuation_table_five['Order'].value_counts(normalize=True)
fraction_ten_order = total_fluctuation_table_ten['Order'].value_counts(normalize=True)
fraction_twenty_order = total_fluctuation_table_twenty['Order'].value_counts(normalize=True)
fraction_fifty_order = total_fluctuation_table_fifty['Order'].value_counts(normalize=True)
fraction_hundred_order = total_fluctuation_table_hundred['Order'].value_counts(normalize=True)
fraction_all_order = total_fluctuation_table['Order'].value_counts(normalize=True)

# Get the fraction of substitutions in each family at decreasing TSS

fraction_five_family = total_fluctuation_table_five['Family'].value_counts(normalize=True)
fraction_ten_family = total_fluctuation_table_ten['Family'].value_counts(normalize=True)
fraction_twenty_family = total_fluctuation_table_twenty['Family'].value_counts(normalize=True)
fraction_fifty_family = total_fluctuation_table_fifty['Family'].value_counts(normalize=True)
fraction_hundred_family = total_fluctuation_table_hundred['Family'].value_counts(normalize=True)
fraction_all_family = total_fluctuation_table['Family'].value_counts(normalize=True)

# Prepare dataframe with fraction in each order versus TSS limit

fraction_across_order = pd.DataFrame({'All': fraction_all_order, '100': fraction_hundred_order, '50':fraction_fifty_order, '20':fraction_twenty_order, '10':fraction_ten_order, '5':fraction_five_order})
fraction_across_order = fraction_across_order.fillna(0)
fraction_across_order = fraction_across_order.drop('MIXED')
fraction_across_order = fraction_across_order.drop('Squamata')

# Prepare dataframe with fraction in each family versus TSS limit

fraction_across_family = pd.DataFrame({'All': fraction_all_family, '100': fraction_hundred_family, '50':fraction_fifty_family, '20':fraction_twenty_family, '10':fraction_ten_family, '5':fraction_five_family})
fraction_across_family = fraction_across_family.fillna(0)
fraction_across_family = fraction_across_family.drop('MIXED')
fraction_across_family = fraction_across_family.drop('Dactyloidae')

# Calculate log2 fold change at each TSS cutoff versus the fraction distributed across orders at all TSS values

fraction_across_order['Log2_Fold_5_vs_All'] = np.log2(fraction_across_order['5']/fraction_across_order['All'])
fraction_across_order['Log2_Fold_10_vs_All'] = np.log2(fraction_across_order['10']/fraction_across_order['All'])
fraction_across_order['Log2_Fold_20_vs_All'] = np.log2(fraction_across_order['20']/fraction_across_order['All'])
fraction_across_order['Log2_Fold_50_vs_All'] = np.log2(fraction_across_order['50']/fraction_across_order['All'])
fraction_across_order['Log2_Fold_100_vs_All'] = np.log2(fraction_across_order['100']/fraction_across_order['All'])

# Calculate log2 fold change at each TSS cutoff versus the fraction distributed across families at all TSS values

fraction_across_family['Log2_Fold_5_vs_All'] = np.log2(fraction_across_family['5']/fraction_across_family['All'])
fraction_across_family['Log2_Fold_10_vs_All'] = np.log2(fraction_across_family['10']/fraction_across_family['All'])
fraction_across_family['Log2_Fold_20_vs_All'] = np.log2(fraction_across_family['20']/fraction_across_family['All'])
fraction_across_family['Log2_Fold_50_vs_All'] = np.log2(fraction_across_family['50']/fraction_across_family['All'])
fraction_across_family['Log2_Fold_100_vs_All'] = np.log2(fraction_across_family['100']/fraction_across_family['All'])

# Save results in CSV files

fraction_across_order.to_csv('Fraction_Across_Order_Log2_Fold_All_100_50_20_10_5.csv')
fraction_across_family.to_csv('Fraction_Across_Family_Log2_Fold_All_100_50_20_10_5.csv')

