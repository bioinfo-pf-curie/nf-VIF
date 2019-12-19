# nf-VIF: A Nextflow-based Virus Insertion Finder

**Institut Curie - Bioinformatics Core Facility**

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)
[![MultiQC](https://img.shields.io/badge/MultiQC-1.6-blue.svg)](https://multiqc.info/)
[![Install with](https://anaconda.org/anaconda/conda-build/badges/installer/conda.svg)](https://conda.anaconda.org/anaconda)
[![Singularity Container available](https://img.shields.io/badge/singularity-available-7E4C74.svg)](https://singularity.lbl.gov/)

### Introduction

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple computing infrastructures in a very portable manner.
It comes with conda / singularity containers making installation trivial and results highly reproducible, and can be run on a single laptop as well as on a cluster.

The current workflow is based on the nf-core best practice. See the nf-core project from details on [guidelines](https://nf-co.re/).

### Pipline summary

This pipeline was designed to process Illumina sequencing data from the HPV capture protocol.
Briefly, it allows to detect and genotype the HPV strain(s) available in the samples, and to precisely map the insertion sites on the Human genome.

1. Reads cleaning and qality controls ([TrimGalore!](https://github.com/FelixKrueger/TrimGalore), [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. HPV Genotyping (([Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)))
2. Local alignment on detected HPV strain(s) ([Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml))
3. Detection of putative HPV breakpoints using soft-clipped reads
4. Soft-clipped reads alignment on Human genome reference ([BLAT](https://genome.ucsc.edu/FAQ/FAQblat.html))
5. Detection of insertion loci and filtering of the results
6. Presentation of results in a dynamic report ([MultiQC](http://multiqc.info/))

### Quick help

```bash
N E X T F L O W  ~  version 19.04.0
Launching `main.nf` [backstabbing_roentgen] - revision: 93bf83bb3b
HPV v1.0dev
=======================================================

Usage:

nextflow run main.nf --reads '*_R{1,2}.fastq.gz' --genome 'hg19'
nextflow run main.nf --samplePlan sample_plan --genome 'hg19'

Mandatory arguments:
  --reads                       Path to input data (must be surrounded with quotes)
  --samplePlan                  Path to sample plan file if '--reads' is not specified
  --genome                      Name of Human genome reference
  -profile                      Configuration profile to use. Can use multiple (comma separated)
                                Available: conda, singularityPath, cluster, test and more.

Options:
  --singleEnd                   Specifies that the input is single end reads

Genome References:              If not specified in the configuration file or you wish to overwrite any of the references.
  --genome                      Name of iGenomes reference
  --bwt2_index                  Path to Bowtie2 index
  --fasta                       Path to Fasta reference (.fasta)
  --blatdb                      Path to BLAT database (.2bit)

HPV References:
  --fasta_hpv                   Path to Fasta HPV reference (.fasta)
  --bwt2_index_hpv              Path to Bowtie2 index for all HPV strains
  --bwt2_index_hpv_split        Path to Bowtie2 index per HPV strain

Other options:
  --nb_geno                     Number of HPV genotype to consider. Default: 3
  --split_report                Generate one report per sample
  --skip_trimming               Skip trimming step
  --skip_fastqc                 Skip quality controls on sequencing reads
  --skip_blat                   Skip Human mapping with Blat
  --skip_multiqc                Skip report
  --outdir                      The output directory where the results will be saved
  --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
  -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
```

### Quick run

The pipeline can be run on any infrastructure from a list of input files or from a sample plan as follow.
Note that by default, all tools are expected to be available from your `PATH`. See the full [`documentation`]('docs/README.md') for details and containers usage.

#### Run the pipeline on a test dataset
See the conf/test.conf to set your test dataset.

```
nextflow run main.nf -profile test

```

#### Run the pipeline from a sample plan

```
nextflow run main.nf --samplePlan MY_SAMPLE_PLAN --genome 'hg19' --outdir MY_OUTPUT_DIR

```

#### Run the pipeline on a cluster

```
echo "nextflow run main.nf --reads '*.R{1,2}.fastq.gz' --genome 'hg19' --outdir MY_OUTPUT_DIR -profile singularityPath,cluster" | qsub -N illumina-hpv

```

### Documentation

1. [Installation](docs/installation.md)
2. [Reference genomes](docs/configuration/reference_genomes.md)  
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)


#### Credits

This pipeline has been set up and written by the sequencing facility, the genetic service and the bioinformatics platform of the Institut Curie (M. Deloger, S. Lameiras, S. Baulande, N. Servant)

#### Contacts

For any question, bug or suggestion, please, contact the bioinformatics core facility

