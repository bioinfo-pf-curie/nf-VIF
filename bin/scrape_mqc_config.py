#!/usr/bin/env python

# Author(s): Nicolas Servant
# Contact: nicolas.servant@curie.fr

# Copyright Institut Curie 2019
# This software is a computer program whose purpose is to analyze high-throughput sequencing data.
# You can use, modify and/ or redistribute the software under the terms of license (see the LICENSE file for more details).
# The software is distributed in the hope that it will be useful, but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND.
# Users are therefore encouraged to test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data.
# The fact that you are presently reading this means that you have had knowledge of the license and that you accept its terms.

""" 
Scripts to create MultiQC config on-the-fly
"""
import argparse
import os

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("filename")
    args = parser.parse_args()

    # Dump to multiQC config
    config_cov='''\
    hpv_cov_{genotype}:
      file_format: 'tsv'
      section_name: '[{genotype}] - Coverage'
      description: 'from HPV local aligment.'
      plot_type: 'linegraph'
      pconfig:
         namespace: 'hpvcov_{genotype}'
         title: '[{genotype}] Alignment Coverage'
         id: 'hpvcov_{genotype}'
         xlab: '{genotype} genomic sequence'
         ylab: 'reads coverage (CPM)'
    '''

    config_bkp='''\
    hpv_bkp_{genotype}:
      file_format: 'tsv'
      section_name: '[{genotype}] - Breakpoints'
      description: 'from HPV local aligment.'
      plot_type: 'linegraph'
      pconfig:
         namespace: 'hpvbkp_{genotype}'
         title: '[{genotype}] - Breakpoints Position'
         id: 'hpvbkp_{genotype}'
         xlab: '{genotype} genomic sequence'
         ylab: 'reads number'
    '''

    sp_cov='''\
    hpv_cov_{genotype}:
      fn: '*{genotype}_covmatrix.mqc'
    '''
 
    sp_bkp='''\
    hpv_bkp_{genotype}:
      fn: '*{genotype}*bkp.mqc'
    '''

    ## load data
    genolist = []
    if os.stat(args.filename).st_size > 0:
        fh = open(args.filename)
        while True:
            line = fh.readline().rstrip()
            if not line:
                break
            else:
                genolist.append(line)
        fh.close()

    ## create custom content
    print ("custom_data:")
    for g in genolist:
        print (config_cov.format(genotype=g))
        print (config_bkp.format(genotype=g))
              
    print ("sp:")
    for g in genolist:
        print (sp_cov.format(genotype=g))
        print (sp_bkp.format(genotype=g))

    print ("report_section_order:")
    spos = 500
    for g in genolist:
        print ("   hpv_cov_" + g + ":")
        print ("      order: -" + str(spos))
        print ("   hpv_bkp_" + g + ":")
        print ("      order: -" + str(spos+1))
        spos += 2
       

