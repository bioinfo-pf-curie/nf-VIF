#!/usr/bin/env python

# Author(s): Marc Deloger and Dimitri Desvillechabrol
# Contact: nicolas.servant@curie.fr

# Copyright Institut Curie 2019
# This software is a computer program whose purpose is to analyze high-throughput sequencing data.
# You can use, modify and/ or redistribute the software under the terms of license (see the LICENSE file for more details).
# The software is distributed in the hope that it will be useful, but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND.
# Users are therefore encouraged to test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data.
# The fact that you are presently reading this means that you have had knowledge of the license & that you accept its terms.

import os
import sys
import re
import pandas as pd
import numpy as np

filename = sys.argv[1]
base = os.path.basename(filename)
base = re.sub(r'\.tsv$', '', base)

table = open('./' + base + '_table.csv', 'w')
table_filtered = open('./' + base + '_table_filtered.csv', 'w')
try:
    data = pd.read_csv(filename, sep="\t")    

    # Set columns name
    data.columns = ['match', 'mismatch', 'repmatch', 'Ns', 'Qgapcount',
                    'Qgapbases', 'Tgapcount', 'Tgapbases', 'strand', 'Qname',
                    'Qsize', 'Qstart', 'Qend', 'Tname', 'Tsize', 'Tstart', 'Tend',
                    'blockcount', 'blockSizes', 'qStarts', 'tStarts']    

    data['sample'] = os.path.splitext(base)[0]    

    # Split query name values
    new = data['Qname'].str.split('|', expand=True)
    data['read'] = new[0]
    data['genotype'] = new[1]
    data['feature'] = new[2]
    data['position'] = new[3]
    data.drop(columns=['Qname'], inplace=True)    

    # Add column with mapped query's proportion
    data['prop'] = (data['Qend'] - data['Qstart']) / data['Qsize']    

    # Switch/Case to determine the insertion chromosome position
    case1 = (data['feature'] == 'left') & (data['strand'] == '-')
    case2 = (data['feature'] == 'left') & (data['strand'] == '+')
    case3 = (data['feature'] == 'right') & (data['strand'] == '+')
    case4 = (data['feature'] == 'right') & (data['strand'] == '-')
    vals = [data['Tstart'], data['Tend'], data['Tstart'], data['Tend']]
    data['chr_position'] = np.select([case1, case2, case3, case4],
                                     vals, default=np.nan)
    data['chr_position'] = data['chr_position'].astype(int)    
    

    # Filter rows having 'Qgapbases <= 5 & Tgapbases <= 5 & prop >0.9'
    data_filtered = data.query('Qgapbases <= 5 & Tgapbases <= 5 & prop > 0.9')    

    # Rename column "Tname" in "chr" and "Qsize" in "size"
    data_filtered_renamed = data_filtered.rename(
        columns={"Tname": "chr", "Qsize": "size"}
    )    

    # Create a key that will be used to group rows and apply a function on them
    data_filtered_renamed['match_key'] = (
        data_filtered_renamed['read'] + "|"
        + data_filtered_renamed['strand'] + "|"
        + data_filtered_renamed['genotype'] + "|"
        + data_filtered_renamed['feature'] + "|"
        + data_filtered_renamed['position'].astype(str) + "|"
        + data_filtered_renamed['chr'] + "|"
        + data_filtered_renamed['chr_position'].astype(str) + "|"
        + data_filtered_renamed['size'].astype(str)
    )    

    # Take the max value of the column 'match' for each group of rows having
    # the same key, and replace all the match values of the group by the max value
    # get max for each  match key
    max_match_by_key = data_filtered_renamed['match'].groupby(
        data_filtered_renamed['match_key']
    ).max()
    data_filtered_renamed = data_filtered_renamed.set_index('match_key')
    data_filtered_renamed.update(max_match_by_key)    

    cleaned = data_filtered_renamed
    cleaned['score'] = 0
    # Create a key that will be used to group rows
    cleaned['left_and_right_key'] = cleaned['genotype'] + "|" + cleaned['chr']    

    # Set the score value to 1 for each group of rows having the same key, if the 
    # 'feature' column of the group contains 2 different values 'left' AND 'right'
    # (let score to 0 if all rows of each group have a unique 'feature' value
    # 'left' OR 'right')
    # get for each feature if they have 'right' and 'left'
    feature_nunique = cleaned['feature'].groupby(
        cleaned['left_and_right_key']
    ).nunique()
    key_have_both = feature_nunique.loc[feature_nunique == 2].index
    cleaned.loc[cleaned['left_and_right_key'].isin(key_have_both), 'score'] = 1    

    # Create a new key to group rows
    cleaned['enoughly_covered_key'] = (
        cleaned['genotype'] + "|"
        + cleaned['feature'] + "|"
        + cleaned['position'].astype(str) + "|"
        + cleaned['chr'] + "|"
        + cleaned['chr_position'].astype(str)
    )
    # Set the score value to 3 for each group of rows having the same key, if the
    # 'score' column is still at 0 AND the group contains at least 2 rows
    # (aka reads)
    count_row = cleaned['enoughly_covered_key'].value_counts()
    multi_row_cleaned = count_row.loc[count_row > 1].index
    only_multi = cleaned.loc[cleaned['enoughly_covered_key'].isin(multi_row_cleaned)]    

    # Set score to 3 if all score is 1
    zero_score = only_multi['score'].groupby(
        only_multi['enoughly_covered_key']
    ).sum()
    zero_index = zero_score.loc[zero_score == 0].index
    cleaned.loc[cleaned['enoughly_covered_key'].isin(zero_index), 'score'] = 3    

    # Set score to 4 if all scores are 1 so the min is 1
    min_score = only_multi['score'].groupby(
        only_multi['enoughly_covered_key']
    ).min()
    one_index = min_score.loc[min_score == 1].index
    cleaned.loc[cleaned['enoughly_covered_key'].isin(one_index), 'score'] = 4    

    # Set the score value to 4 for each group of rows having the same key, if the
    # 'score' column was previously set to 1 AND the group contains at least
    # 2 rows    

    # Create a key that will be used to group rows		
    cleaned['uniquely_mapped_key'] = (
        cleaned['read'] + "|" + cleaned['genotype'] + "|" + cleaned['feature']
    )    

    # add 6 to score to uniquely_mapped_key
    count_row = cleaned['uniquely_mapped_key'].value_counts()
    one_row_index = count_row.loc[count_row == 1].index
    cleaned.loc[cleaned['uniquely_mapped_key'].isin(one_row_index), 'score'] += 6    
    

    table_to_display = cleaned.loc[:, [
        'sample', 'genotype', 'feature', 'score', 'position', 'chr',
        'chr_position', 'match'
    ]]
    # add count column
    table_to_display['counts'] = 0
    table_to_display['key'] = (
        table_to_display['sample'] + "|"
        + table_to_display['genotype'] + "|"
        + table_to_display['feature'] + "|"
        + table_to_display['position'].astype(str) + "|"
        + table_to_display['chr'] + "|"
        + table_to_display['chr_position'].astype(str)
    )    

    # set max score and max match for each group
    table_to_display = table_to_display.set_index('key', drop=False)
    for col in ('score', 'match'):
        max_match_by_key = table_to_display[col].groupby(
            table_to_display['key']
        ).max()
        table_to_display.update(max_match_by_key)    

    # set number of row for each groupby
    count_row = table_to_display['key'].value_counts()
    count_row.name = 'counts'
    table_to_display.update(count_row)
    table_to_display = table_to_display.drop_duplicates()
    table_to_display = table_to_display.drop(columns=['key'])    

    table_to_display.to_csv(table, index=False, float_format='%.0f')

    # filter table to reduce noise
    table_to_display = table_to_display[
        (table_to_display['score'] >= 4)
        & (table_to_display['match'] >= 10)
        & (table_to_display['counts'] >= 2)
    ]    

    table_to_display.sort_values(by=['score', 'counts'], ascending=False, inplace=True)

    # add index for multiQC
    table_to_display = table_to_display.reset_index()
    table_to_display['multiqc_index'] = (
        table_to_display['sample'] + '|' + table_to_display.index.astype(str)
    )    

    table_to_display = table_to_display.reindex(columns=[
        'multiqc_index', 'sample', 'genotype', 'feature', 'score', 'position',
        'chr', 'chr_position', 'match', 'counts'
    ])    

    table_to_display.to_csv(table_filtered, index=False, float_format='%.0f')

except:
    emptydf=pd.DataFrame(columns=['multiqc_index','sample','genotype','feature','score','position','chr','chr_position','match','counts'])
    emptydf.to_csv(path_or_buf=table_filtered,index=False)
    emptydf.to_csv(path_or_buf=table,index=False)
