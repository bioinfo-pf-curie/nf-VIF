#!/usr/bin/env python

# Author(s): Tina Alaeitabar
# Contact: tina.alaeitabar@curie.fr

# Copyright Institut Curie 2020
# This software is a computer program whose purpose is to analyze high-throughput sequencing data.
# You can use, modify and/ or redistribute the software under the terms of license (see the LICENSE file for more details).
# The software is distributed in the hope that it will be useful, but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND.
# Users are therefore encouraged to test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data.
# The fact that you are presently reading this means that you have had knowledge of the license and that you accept its terms.

""" 
Scripts to find overlap between hits and braekpoints region.
"""
#import pandas as pd
#from operator import itemgetter
#import numpy as np
#import sys
#import argparse
#import os

import os
import sys
import re
import pandas as pd
import numpy as np
import argparse


def parsBlatFile(hits_filename):

   try:
	data = pd.read_csv(hits_filename, sep="\t", header=None)
	data.drop_duplicates(inplace=True)
	data.reset_index(drop=True, inplace=True)

	### Set columns name
	data.columns = ['match', 'mismatch', 'repmatch', 'Ns', 'Qgapcount',
                    'Qgapbases', 'Tgapcount', 'Tgapbases', 'strand', 'Qname',
                    'Qsize', 'Qstart', 'Qend', 'Tname', 'Tsize', 'Tstart', 'Tend',
                    'blockcount', 'blockSizes', 'qStarts', 'tStarts']

	###  Split query name values
	new = data['Qname'].str.split('|', expand=True)
	data['read'] = new[0]
	data['genotype'] = new[1]
	data['feature'] = new[2]
	data['position'] = new[3]
	data.drop(columns=['Qname'], inplace=True)

	### Switch/Case to determine the insertion chromosome position
	case1 = (data['feature'] == 'left') & (data['strand'] == '-')
	case2 = (data['feature'] == 'left') & (data['strand'] == '+')
	case3 = (data['feature'] == 'right') & (data['strand'] == '+')
	case4 = (data['feature'] == 'right') & (data['strand'] == '-')
	vals = [data['Tstart'], data['Tend'], data['Tstart'], data['Tend']]
	data['chr_position'] = np.select([case1, case2, case3, case4],
                                  vals, default=np.nan)
	data['chr_position'] = data['chr_position'].astype(int)
	data['chr_position'] = data.apply(lambda x: "{:,}".format(x['chr_position']), axis=1)

	data['prop'] = (data['Qend'] - data['Qstart']) / data['Qsize']
	data['hg19_orientation'] = np.where(data['strand'] == '+', "Forward","Reverse")
       
         ### Filter rows having 'Qgapbases <= 5 & Tgapbases <= 5 & prop >0.9'
	data_filtered = data.query('Qgapbases <= 5 & Tgapbases <= 5 & prop > 0.9')
	data_filtered = data_filtered.drop(columns=['genotype','feature','strand','mismatch','repmatch','Ns','Qgapcount','Qgapbases','Tgapcount','Tgapbases','Qstart','Qend','Tsize','blockcount','blockSizes','qStarts','Tend','Tstart','tStarts'])

	data_filtered = data_filtered[['read','Tname','chr_position','hg19_orientation','match','Qsize','prop']]
	data_filtered_grouped = data_filtered.groupby(['read','Tname', 'chr_position','hg19_orientation']).agg({'match':'max'})
	data_filtered_grouped = data_filtered_grouped.reset_index()
	data_filtered_grouped = data_filtered_grouped.rename(columns={'match':'match_max'})
	data_filtered = pd.merge(data_filtered, data_filtered_grouped, how='left', on=['read', 'Tname', 'chr_position','hg19_orientation'])
	data_filtered = data_filtered[data_filtered['match'] == data_filtered['match_max']]
        data_filtered.drop(columns=['match'],inplace=True)
	data_filtered.rename(columns={'match_max': 'match'},inplace=True)
	data_filtered.drop_duplicates(inplace=True)
        del data_filtered_grouped
   except:
            print >> sys.stderr, "Error in blat file parsing!", hits_filename;
   return(data_filtered) 


def getMaxFlanqSeq(BP_filename):
   """select reads with maximum flanking_seq"""

   try:
	df=pd.read_csv(BP_filename,sep = '\t')
	#df.drop_duplicates(inplace=True)
	#df.reset_index(drop=True, inplace=True)
	df['hpv_BP_orientation'] = np.where(df['Label'] == '3prime', "<---","--->")
	df['count'] = df.groupby('Sequence')['Sequence'].transform('count')

	df_grouped = df.groupby(['Sequence']).agg({'Length_Seq':'max'})
	df_grouped = df_grouped.reset_index()
	df_grouped = df_grouped.rename(columns={'Length_Seq':'Length_Max'})
	df = pd.merge(df,df_grouped, how='left', on=['Sequence'])
	df = df.loc[df['Length_Seq'] == df['Length_Max']]
      
	flanking_seq=df.groupby('Sequence')['Flanking_Seq'].nth(0)
	df= df.set_index('Sequence')
	df.update(flanking_seq)

	df.sort_values(by ='count', ascending=False, inplace=True)
	df.reset_index(inplace=True)
	df.rename(columns={'Pos_Ref': 'Hpv_BP_position' ,'Flanking_Seq':'flanking_seq'}, inplace=True)
        del df_grouped

   except:
            print >> sys.stderr, "Error in parsing:", BP_filename;
   return (df)


def merged(df_bp_selected,df_hits):
  
   """merged filter hits & reads with longest flanking_seq"""

   try:
	merged = pd.merge(df_bp_selected,df_hits, how='left', right_on=['read'], left_on = ['Read_Name'],indicator=True)
	merged.drop_duplicates(inplace=True)
	merged.drop(columns=['Read_Name','Coordinate_Fasta','Length_Seq','Breakpoint','Qsize','SoftClip_Position','read','BraekPoint_Position','Length_Max','CigarString'], inplace=True)
	merged_both=merged.loc[merged['_merge'] =='both'].copy()
	merged_left=merged.loc[merged['_merge'] =='left_only'].copy()
	del[merged]

	### Find hits & multi_hits
	merged_both['match'] = merged_both.groupby(['Sequence'])['match'].transform(max)
	merged_both['chr_position'] = merged_both.groupby(['Sequence'])['chr_position'].transform(max)
	merged_both.drop_duplicates(inplace=True)

	merged_both['chromosom_count'] = merged_both.groupby('Sequence')['Tname'].transform('count')

	hg19_orientation=merged_both.groupby(['Sequence'], as_index = True).agg({'hg19_orientation': '|'.join})
	merged_both = merged_both.set_index('Sequence')
	merged_both.update(hg19_orientation)

	merged_both.reset_index(inplace=True)
	Tname=merged_both.groupby(['Sequence'], as_index = True).agg({'Tname': '|'.join})
	merged_both = merged_both.set_index('Sequence')
	merged_both.update(Tname)


	merged_both['hg19_BP_chromosom'] = np.where(merged_both['chromosom_count'] == 1, 'uniq','Multi_hit')
	merged_both['hg19_BP_chromosom'] = np.where(merged_both['hg19_BP_chromosom'] == 'uniq', merged_both['Tname'],'Multi_hit')

	merged_both['hg19_orientation'] = np.where(merged_both['chromosom_count'] == 1, merged_both['hg19_orientation'],'-')
	merged_both['hg19_orientation'] = np.where(merged_both['hg19_orientation'] == 'forward', '--->','<---')
	merged_both['hg19_BP_position'] = np.where(merged_both['chromosom_count'] == 1, merged_both['chr_position'],'-')

	merged_both_order = merged_both[['count','Hpv_BP_position','hpv_BP_orientation','hg19_BP_chromosom','hg19_BP_position','hg19_orientation','match','flanking_seq']]
        del merged_both

	### Find no_hits
	merged_left['match'] = 0
	merged_left['hg19_BP_chromosom'] = 'No_hits'
	merged_left['hg19_BP_position'] = '-'
	merged_left['hg19_orientation'] = '-'
	merged_left['Tname'] = '-'

	merged_left = merged_left.set_index('Sequence')
	merged_left_order = merged_left[['count','Hpv_BP_position','hpv_BP_orientation','hg19_BP_chromosom','hg19_BP_position','hg19_orientation','match','flanking_seq']]
        del merged_left

   except:
        print >> sys.stderr, "Error in merge:",;
   return (merged_both_order, merged_left_order)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--blat_file", help="Hits file.", default="")
    parser.add_argument("-b", "--bp_file", help="Breaking Point file", default="")
    parser.add_argument("-o", "--outputDir", help="Output directory", default="./")
    args = parser.parse_args()

    if args.blat_file == "":
       print ('Blat file is required.')
       exit(0) 
   
    if args.bp_file == "":
       print ('Breakpoint are required.')
       exit(0) 

    try:
        hits_filename=args.blat_file
        bp_filename=args.bp_file
        base = os.path.basename(hits_filename)
        base = re.sub(r'\.tsv$', '', base)

	df_hits=parsBlatFile(hits_filename)
        df_bp_selected=getMaxFlanqSeq(bp_filename)
        (merged_both_order, merged_left_order)=merged(df_bp_selected,df_hits)
        del df_bp_selected
        del df_hits

        ### concatenate hits & no_hit
        concat=pd.concat([merged_both_order,merged_left_order])
        concat.sort_values(by ='count', ascending=False, inplace=True)
        concat.drop_duplicates(inplace=True)

        concat['Hpv_BP_position'] = concat.apply(lambda x: "{:,}".format(x['Hpv_BP_position']), axis=1)
        concat.insert(0,'sample',os.path.splitext(base)[0])
        concat.rename(columns={'count':'Count','Hpv_BP_position': 'HPV_BP_position','hpv_BP_orientation': 'HPV_BP_orientation', 'match': 'Match', 'Flanking_Seq':'flanking_seq'}, inplace=True)
        
        # add index for multiQC
        concat.reset_index(inplace=True)
        concat.insert(0,'multiqc_index', concat['sample'] + '|' + concat.index.astype(str))
        concat.to_csv('./' + base + '_table.csv', index=False, sep = "\t", float_format='%.0f')

        # filter rows to reduce noise
        concat_filtered=concat[concat['Count']>1]
        concat_filtered.to_csv('./' + base + '_table_filtered.csv', index=False, sep = "\t", float_format='%.0f')

    except:
        print >> sys.stderr, "Error in main:\n";
        emptydf=pd.DataFrame(columns=['Sequence','sample','Count','HPV_BP_position','HPV_BP_orientation','hg19_BP_chromosom','hg19_BP_position','hg19_orientation','Match''flanking_seq'])
        emptydf.to_csv(path_or_buf=table_filtered,index=False)
        emptydf.to_csv(path_or_buf=table,index=False)

