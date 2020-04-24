#!/usr/bin/env python

# Author(s): Marc Deloger, Dimitri Desvillechabrol, Tina Alaeitabar, Nicolas Servant
# Contact: nicolas.servant@curie.fr

# Copyright Institut Curie 2019-2020
# This software is a computer program whose purpose is to analyze high-throughput sequencing data.
# You can use, modify and/ or redistribute the software under the terms of license (see the LICENSE file for more details).
# The software is distributed in the hope that it will be useful, but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND.
# Users are therefore encouraged to test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data.
# The fact that you are presently reading this means that you have had knowledge of the license & that you accept its terms.

import argparse
import os
import sys
import re
import pandas as pd
import numpy as np

pd.set_option('display.max_rows', None, 'display.max_columns', None, 'display.max_colwidth', 200)


def loadBlat(filename, minProp=0.8):
    """
    Load and filter BLAT hits file
    """
    data = pd.read_csv(filename, sep="\t", header=None)    
    base = os.path.basename(filename)
        
    # Set columns name
    data.columns = ['match', 'mismatch', 'repmatch', 'Ns', 'Qgapcount',
                    'Qgapbases', 'Tgapcount', 'Tgapbases', 'strand', 'Qname',
                    'Qsize', 'Qstart', 'Qend', 'Tname', 'Tsize', 'Tstart', 'Tend',
                    'blockcount', 'blockSizes', 'qStarts', 'tStarts']    

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
    end_vals = [data['Tend'], data['Tstart'], data['Tend'], data['Tstart']]
    data['chr_position'] = np.select([case1, case2, case3, case4], vals, default=np.nan)
    data['chr_position'] = data['chr_position'].astype(int)    
    data['end'] = np.select([case1, case2, case3, case4], end_vals, default=np.nan)
    data['end'] = data['end'].astype(int)
    
    ## update colums
    data = data.rename(columns={"Tname": "chr", "Qsize": "size"})
    data['human_orientation'] = np.where(data['strand'] == '+', "forward","reverse")
    data['virus_strand'] = np.where(data['feature'] == '3prime', "<---","--->")

    # Filter rows having 'Qgapbases <= 5 & Tgapbases <= 5 & prop >0.9'
    data = data.query('Qgapbases <= 5 & Tgapbases <= 5 & prop > @minProp')

    ## Remove unused columns
    data = data.drop(columns=['mismatch','repmatch','Ns','Qgapcount','Qgapbases','Tgapcount','Tgapbases','prop',
                              'Qstart','Qend','Tsize','blockcount','blockSizes','qStarts','Tend','Tstart','tStarts'])
    
    return(data)


def groupBlat(data_grouped):
    """
    Combine rows in the BLAT table
    """ 
    # Filter duplicates rows to deal with R1==R2 reads
    data_grouped.drop_duplicates(inplace=True)

    #Create a key that will be used to group rows and apply a function on them
    data_grouped['match_key'] = (
        data_grouped['read'] + "|" 
        + data_grouped['strand'] + "|"
        + data_grouped['genotype'] + "|"
        + data_grouped['feature'] + "|"
        + data_grouped['position'].astype(str) + "|"
        + data_grouped['chr'] + "|"
        + data_grouped['chr_position'].astype(str)
        + data_grouped['size'].astype(str)
    )
        
    # Take the max value of the column 'match' for each group of rows having
    # the same key, and replace all the match values of the group by the max value
    # get max for each match key
    max_match_by_key = data_grouped['match'].groupby(
        data_grouped['match_key']
    ).max()
    data_grouped = data_grouped.set_index('match_key')
    data_grouped.update(max_match_by_key)

    # Create a key to define a single breakpoint
    data_grouped['bkp_key'] = (
        data_grouped['genotype'] + "|"
        + data_grouped['feature'] + "|"
        + data_grouped['position'].astype(str) + "|"
        + data_grouped['chr'] + "|"
        + data_grouped['strand'] + "|"
        + data_grouped['chr_position'].astype(str)
    )
    
    # For a given breakpoint, we can have multiple read ends.
    # Report the max one
    ends_by_key = data_grouped['end'].groupby(
        data_grouped['bkp_key']
    ).max()
    data_grouped = data_grouped.set_index('bkp_key', drop=False)
    data_grouped.update(ends_by_key)
      
    ## set max match for each breakpoint
    max_match_by_key = data_grouped['match'].groupby(
        data_grouped['bkp_key']
    ).max()
    data_grouped.update(max_match_by_key)

    ## Add a count column with the number of reads for a given breakpoint
    ## Use the number of unique read name
    data_grouped['countHQ'] = 0
    data_grouped['countHQ'] = data_grouped['read'].groupby(
        data_grouped['bkp_key']
    ).nunique()
     
    # remove key columns and index
    data_grouped = data_grouped.drop(columns=['bkp_key'])
    data_grouped = data_grouped.reset_index(drop=True)

    return(data_grouped)
    

def getMaxFlanqSeq(BP_filename):
    """
    Select reads with maximum flanking_seq
    Read_Name Sequence Pos_Ref CigarString Label SoftClip_Position BreakPoint_Position Breakpoint Coordinate_Fasta Length_Seq Flanking_Seq
    """
    df=pd.read_csv(BP_filename,sep = '\t')

    # Filter duplicates rows to deal with R1==R2 reads
    df.drop_duplicates(inplace=True)

    ## Counts unique read names with the same motif
    df['count']=df.groupby('Sequence')['Read_Name'].transform('nunique') 
    df['virus_bkp_orientation'] = np.where(df['Label'] == '3prime', "<---","--->")
 
    ## Add the Max_Flanking_Seq for each cases
    df = df.set_index('Sequence')
    df['Length_Max'] = df.groupby(['Sequence'])['Length_Seq'].transform('max')
    df_max = df.loc[df['Length_Seq'] == df['Length_Max']]
    df_max = df_max[['Flanking_Seq']]

    ## Select randomly one max sequence if multiple available
    fseq = df_max.groupby(['Sequence'])['Flanking_Seq'].nth(0)
    df_max.update(fseq)
    df_max = df_max.rename(columns={'Flanking_Seq':'max_flanking_seq'})
    df = pd.merge(df,df_max, how='left', on=['Sequence'])
   
    #df_grouped = df.groupby(['Sequence']).agg({'Length_Seq':'max'})
    #df_grouped = df_grouped.reset_index()
    #df_grouped = df_grouped.rename(columns={'Length_Seq':'Length_Max'})
    #df = pd.merge(df,df_grouped, how='left', on=['Sequence'])
    #df = df.loc[df['Length_Seq'] == df['Length_Max']]
    #flanking_seq=df.groupby('Sequence')['Flanking_Seq'].nth(0)
    #df= df.set_index('Sequence')
    #df.update(flanking_seq)

    df.sort_values(by ='count', ascending=False, inplace=True)
    df.reset_index(inplace=True)
    df.rename(columns={'Pos_Ref': 'virus_bkp_pos'}, inplace=True)
    del df_max

    return (df)
    
    
def merged(df_bp_selected, df_hits, hideMultihits=False):
    """
    Merged BLAT hits & flanking_seq information
    The main key here is the motif around the breakpoint
    """

    ## Need to take into account the virus position to merge
    df_bp_selected['mkey'] = (df_bp_selected['Read_Name'] + "|" + df_bp_selected['virus_bkp_pos'].astype(str))
    df_hits['mkey'] = (df_hits['read'] + "|" + df_hits['position'].astype(str))
    merged = pd.merge(df_bp_selected, df_hits, how='left', right_on=['mkey'], left_on = ['mkey'], indicator=True)
    #merged = pd.merge(df_bp_selected,df_hits, how='left', right_on=['read'], left_on = ['Read_Name'],indicator=True)

    merged.drop(columns=['Read_Name','Coordinate_Fasta','Length_Seq','Length_Max','Breakpoint','size',
                         'Flanking_Seq','position','Label','strand',
                         'SoftClip_Position','read','BreakPoint_Position','CigarString'], inplace=True)
    merged_both=merged.loc[merged['_merge'] =='both'].copy()
    merged_left=merged.loc[merged['_merge'] =='left_only'].copy()
    del[merged]
    merged_both.drop_duplicates(inplace=True)
    merged_left.drop_duplicates(inplace=True)

    ## Count number of mapping hits using both chromosome and position
    merged_both['chr_key'] = (
        merged_both['chr'] + "|"
        + merged_both['chr_position'].astype(str)
    )
    merged_both['chr_pos_count'] = merged_both.groupby(['Sequence'])['chr_key'].transform('nunique')

    ## In case of clustered integration, get range of insertion
    merged_both['chr_pos_max'] =  merged_both.groupby(['Sequence'])['chr_position'].transform('max')
    merged_both['chr_pos_min'] = merged_both.groupby(['Sequence'])['chr_position'].transform('min')
    merged_both['chr_pos_range'] = '[' + merged_both['chr_pos_min'].astype(int).map(str) + '-' + merged_both['chr_pos_max'].astype(int).map(str) + ']'
    merged_both['end_max'] = merged_both.groupby(['Sequence'])['end'].transform('max')
    merged_both['end_min'] = merged_both.groupby(['Sequence'])['end'].transform('min') 
    merged_both['end_range'] = '[' + merged_both['end_min'].astype(int).map(str) + '-' + merged_both['end_max'].astype(int).map(str) + ']'

    merged_both['chr_count'] = merged_both.groupby(['Sequence'])['chr'].transform('nunique')
    merged_both = merged_both.drop(columns=['chr_key', 'mkey', 'chr_pos_min', 'chr_pos_max', 'end_min', 'end_max'])

    ## Match and Score value
    merged_both['match'] = merged_both.groupby(['Sequence'])['match'].transform(max)
    merged_both['score'] = merged_both.groupby(['Sequence'])['score'].transform(max)
    merged_both['countHQ'] = merged_both.groupby(['Sequence'])['countHQ'].transform(max) ## not sure ...
    merged_both = merged_both.set_index('Sequence')

    ## Chromosomes
    Tname = merged_both.groupby(['Sequence'], as_index = True)['chr'].apply(lambda x: "%s" % '|'.join(np.unique((x))))
    merged_both.update(Tname)
    
    ## Transform multiHits information to ease table interpretation
    merged_both['chr'] = np.where(merged_both['chr_pos_count'] == 1, merged_both['chr'], 
                                  np.where(merged_both['chr_count'] == 1, merged_both['chr'], 'multiHits'))

    merged_both['human_orientation'] = np.where(merged_both['chr_pos_count'] == 1, merged_both['human_orientation'],'-')
    merged_both['human_orientation'] = np.where(merged_both['human_orientation'] == 'forward', '--->','<---')

    merged_both['chr_position'] = np.where(merged_both['chr_pos_count'] == 1, merged_both['chr_position'].astype(int), 
                                           np.where(merged_both['chr_count'] == 1, merged_both['chr_pos_range'], '-'))
    merged_both['end'] = np.where(merged_both['chr_pos_count'] == 1, merged_both['end'].astype(int),
                                  np.where(merged_both['chr_count'] == 1, merged_both['end_range'], '-'))

    merged_both = merged_both.rename(columns={"chr":"human_bkp_chr", "chr_position":"human_bkp_pos"})

    merged_both.drop_duplicates(inplace=True)

    ## no Hits
    merged_left['match'] = 0
    merged_left['human_bkp_chr'] = 'noHits'
    merged_left['human_bkp_pos'] = ''
    merged_left['human_orientation'] = ''
    merged_left['chr'] = ''
    merged_left = merged_left.set_index('Sequence')

    ## Export
    merged_both_order = merged_both[['genotype', 'feature', 'countHQ', 'count', 'score', 'virus_bkp_pos', 'virus_bkp_orientation', 
                                     'human_bkp_chr', 'human_bkp_pos', 'end', 'human_orientation', 'match', 'max_flanking_seq']]
    
    merged_left_order = merged_left[['genotype', 'feature', 'countHQ', 'count', 'score', 'virus_bkp_pos', 'virus_bkp_orientation',
                                     'human_bkp_chr', 'human_bkp_pos', 'end', 'human_orientation', 'match', 'max_flanking_seq']]

    ## Clean
    del merged_both
    del merged_left

    return(merged_both_order, merged_left_order)
    

def addScore(df):
    """
    Add score to classify insertion sites
    """
    df['score'] = 0

    ###################
    ## +1 is at least one left and one rigth site are detected on the same chromosome
    df['left_and_right_key'] = df['genotype'] + "|" + df['chr']
    feature_nunique = df['feature'].groupby(
        df['left_and_right_key']
    ).nunique()
    key_have_both = feature_nunique.loc[feature_nunique == 2].index
    df.loc[df['left_and_right_key'].isin(key_have_both), 'score'] = 1    

    ###################
    ## +3 if at least 2 reads support the same breakpoint
    df['bkp_key'] = (
        df['genotype'] + "|"
        + df['feature'] + "|"
        + df['position'].astype(str) + "|"
        + df['chr'] + "|"
        + df['chr_position'].astype(str)
    )
    # Set the score value to 3 for each group of rows having the same key, if the
    # 'score' column is still at 0 AND the group contains at least 2 reads
    count_row = df['bkp_key'].value_counts()
    multi_row_cleaned = count_row.loc[count_row > 1].index
    only_multi = df.loc[df['bkp_key'].isin(multi_row_cleaned)]

    # Set score to 3 if all score is 1
    zero_score = only_multi['score'].groupby(
        only_multi['bkp_key']
    ).sum()
    zero_index = zero_score.loc[zero_score == 0].index
    df.loc[df['bkp_key'].isin(zero_index), 'score'] = 3    

    # Set score to 4 if all scores are 1 so the min is 1
    min_score = only_multi['score'].groupby(
        only_multi['bkp_key']
    ).min()
    one_index = min_score.loc[min_score == 1].index
    df.loc[df['bkp_key'].isin(one_index), 'score'] = 4    

    #############
    ## +6 is the read align uniquely
    df['uniquely_mapped_key'] = (
        df['read'] + "|" + df['genotype'] + "|" + df['feature']
    )    

    # add 6 to score to uniquely_mapped_key
    count_row = df['uniquely_mapped_key'].value_counts()
    one_row_index = count_row.loc[count_row == 1].index
    df.loc[df['uniquely_mapped_key'].isin(one_row_index), 'score'] += 6    

    ############
    ## Set median score per breakpoint
    med_score_by_key = df['score'].groupby(
        df['bkp_key']
    ).median()
    df = df.set_index('bkp_key', drop=False)
    df.update(med_score_by_key)

    # remove keys
    df = df.drop(columns=['left_and_right_key']) 
    df = df.drop(columns=['uniquely_mapped_key'])
    df = df.drop(columns=['bkp_key'])
    df = df.reset_index(drop=True)

    return (df)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--blat", help="BLAT output file")
    parser.add_argument("-b", "--bkp", help="Breakpoint file with motif and flanking sequence", default=None)
    parser.add_argument("-p", "--minProp", help="Minimum BLAT mapping proportion to keep a read", default=0.8)
    parser.add_argument("-o", "--outputDir", help="Output directory", default="./")
    parser.add_argument("-s", "--sname", help="Sample Name", default="")

    args = parser.parse_args()

    ## I/O
    base = os.path.basename(args.blat)
    base = re.sub(r'\.tsv$', '', base)
    table = open(args.outputDir + "/" + base + '_table.csv', 'w')
    table_filtered = open(args.outputDir + "/" + base + '_table_filtered.csv', 'w')
    if args.bkp:
        tableBkp = open(args.outputDir + "/" + base + '_bkptable.csv', 'w')
        tableBkp_filtered = open(args.outputDir + "/" + base + '_bkptable_filtered.csv', 'w')

    try:
        ##########################################
        ## BLAT results
        df_hits = loadBlat(args.blat, minProp=args.minProp)
        df_grouped = groupBlat(df_hits)
        df_grouped = addScore(df_grouped)
 
        if args.sname:
            df_grouped.insert(0,'sample', args.sname)
        else:
            df_grouped.insert(0,'sample',os.path.splitext(base)[0])

        table_to_display = df_grouped.reindex(columns=[
            'multiqc_index', 'sample', 'genotype', 'feature', 'score', 'position','virus_strand',
            'chr', 'chr_position', 'end', 'strand', 'match', 'countHQ'])
        table_to_display = table_to_display.drop_duplicates()
        table_to_display = table_to_display.reset_index(drop=True)
        table_to_display['multiqc_index'] = (
            table_to_display['sample'] + '|' + table_to_display.index.astype(str)
        )
        table_to_display.sort_values(by=['score', 'countHQ'], ascending=False, inplace=True)
        table_to_display = table_to_display.rename(columns={"position": "virus_bkp_pos", "chr":"human_bkp_chr",
                                                            "chr_position":"human_bkp_pos", "end":"human_end",
                                                            "strand":"human_orientation", "virus_strand":"virus_bkp_orientation"})

        table_to_display.to_csv(path_or_buf=table, index=False, sep =",", float_format='%.0f')

        table_to_display = table_to_display[
            (table_to_display['score'] >= 4)
            & (table_to_display['match'] >= 10)
            & (table_to_display['countHQ'] >= 2)
        ]
        
        if table_to_display.empty:
            table_to_display.append(pd.Series([np.nan]), ignore_index=True)
        table_to_display.to_csv(path_or_buf=table_filtered, index=False, sep =",", float_format='%.0f')

        ########################################
        ## Merge with breakpoint information and flanking sequences
        if args.bkp:
            df_bp_selected = getMaxFlanqSeq(args.bkp)
            (merged_both_order, merged_left_order)=merged(df_bp_selected, df_grouped)
            del df_bp_selected
            del df_hits

            ### concatenate hits & no_hit
            concat=pd.concat([merged_both_order,merged_left_order])
            concat.sort_values(by='count', ascending=False, inplace=True)
            concat.drop_duplicates(inplace=True)

            # add sample name
            if args.sname:
                concat.insert(0,'sample',args.sname) 
            else:
                concat.insert(0,'sample',os.path.splitext(base)[0])
            
            # add index for multiQC
            concat.reset_index(inplace=True)
            concat.insert(0,'multiqc_index', concat['sample'] + '|' + concat.index.astype(str))
            # save
            concat.reset_index()
            concat = concat.rename(columns={"Flanking_Seq": "flanking_seq", "Sequence":"motif", "end":"human_end"}) 
            concat.to_csv(path_or_buf=tableBkp, index=False, sep = ",", float_format='%.0f')
        
            ### Filtered results
            # filter rows to reduce noise
            concat_filtered=concat[(concat['count']>=2) & (concat['countHQ']>=2)]
            concat_filtered=concat_filtered[concat_filtered['human_bkp_chr'] != 'noHits' ]

            ## Add an empty row if no data
            if concat_filtered.empty:
                concat_filtered=concat_filtered.append(pd.Series([np.nan]), ignore_index=True)
            concat_filtered.to_csv(path_or_buf=tableBkp_filtered, index=False, sep = ",", float_format='%.0f')

    except pd.io.common.EmptyDataError:
        sys.stderr.write("Input files is empty and has been skipped.")
        emptydf=pd.DataFrame(columns=['genotype', 'feature', 'countHQ', 'count', 'score', 'virus_bkp_pos', 'virus_bkp_orientation',
                                      'human_bkp_chr', 'human_bkp_pos', 'end', 'human_orientation', 'match', 'Flanking_Seq'])
        emptydf.to_csv(path_or_buf=table_filtered,index=False)
        emptydf.to_csv(path_or_buf=table,index=False)

    except:
        raise

        
