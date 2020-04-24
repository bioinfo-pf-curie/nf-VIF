**********************************
## CHANGES IN VERSION 1.1.1

NEW FEATURES

   o First test to infer the type of sample, EPI/INTEGRATED/NEG

   o Add General Metrics table with additional info

   o Move to Python3

SIGNIFICANT USER-VISIBLE CHANGES

   o Display insertion site table even if no breakpoint is detected

   o Remove non useful Fastqc plots

**********************************
## CHANGES IN VERSION 1.1.0

NEW FEATURES

   o New --minMapq option to provide the minimum mapping quality for the end-to-end mapping steps

   o New --minLen option to provide the minimum length of clipped sequence to consider

   o Update of BLAT summary with flanking sequence info
   
   o Manage overlapping R1/R2 reads reporting the same breakpoint

SIGNIFICANT USER-VISIBLE CHANGES

   o Report genotype available at --FreqGeno % of hpv mapped reads

   o Coverage plot - y axis in reads count when splitReport is specified (CPM otherwise)

   o Report end as max reads position

***********************************
## CHANGES IN VERSION 1.0.1

SIGNIFICANT USER-VISIBLE CHANGES

   o Add BLAT Strand orientation in summary table

   o Change minLen parameter to 10 for extractSoftclipped.py

   o Filter clipped reads based on minLen parameter before breakpoint frequency calculation

BUG FIXES

   o Fix a major bug in nf channels

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
