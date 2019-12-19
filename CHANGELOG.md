***********************************
## CHANGES IN VERSION 1.0.0

NEW FEATURES

  o Add BLAT database in the params.genome

  o Add TrimGalore! for reads trimming

  o Add option to distinguish 3'/5' bkp

  o Add python script to extract soft-clipped reads

  o Add MultiQC report

  o Add mapping on control regions

  o Add Bowtie2 (global and local) mapping on HPV strains

  o First Nextflow version which implements the analysis strategy for HPV analysis from S. Lameiras et al.

BUG FIXES

  o Fix bug when HPV ref is provided in command line (#12)

  o Deal with cases with no detected genotypes
