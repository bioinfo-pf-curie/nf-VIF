**********************************
version-1.1.2

SIGNIFICANT USER-VISIBLE CHANGES

   - Export sorted BAM/BAI files

BUG FIXES

   - Use CPM normalization only when all samples are merged in the same report (#13)

**********************************
version-1.1.1

NEW FEATURES

   - First test to infer the type of sample, EPI/INTEGRATED/NEG
   - Add General Metrics table with additional info
   - Move to Python3

SIGNIFICANT USER-VISIBLE CHANGES

   - Display insertion site table even if no breakpoint is detected
   - Remove non useful Fastqc plots

**********************************
version-1.1.0

NEW FEATURES

   - New --minMapq option to provide the minimum mapping quality for the end-to-end mapping steps
   - New --minLen option to provide the minimum length of clipped sequence to consider
   - Update of BLAT summary with flanking sequence info
   - Manage overlapping R1/R2 reads reporting the same breakpoint

SIGNIFICANT USER-VISIBLE CHANGES

   - Report genotype available at --FreqGeno % of hpv mapped reads
   - Coverage plot - y axis in reads count when splitReport is specified (CPM otherwise)
   - Report end as max reads position

***********************************
version-1.0.1

SIGNIFICANT USER-VISIBLE CHANGES

   - Add BLAT Strand orientation in summary table
   - Change minLen parameter to 10 for extractSoftclipped.py
   - Filter clipped reads based on minLen parameter before breakpoint frequency calculation

BUG FIXES

   - Fix a major bug in nf channels

***********************************
version-1.0.0

NEW FEATURES

  - Add BLAT database in the params.genome
  - Add TrimGalore! for reads trimming
  - Add option to distinguish 3'/5' bkp
  - Add python script to extract soft-clipped reads
  - Add MultiQC report
  - Add mapping on control regions
  - Add Bowtie2 (global and local) mapping on HPV strains
  - First Nextflow version which implements the analysis strategy for HPV analysis from S. Lameiras et al.

BUG FIXES

  - Fix bug when HPV ref is provided in command line (#12)
  - Deal with cases with no detected genotypes
