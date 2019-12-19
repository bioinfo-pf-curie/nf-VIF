# Output

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

## Pipeline overview
The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [Quality Controls](#trimgalore) - reads trimming and quality control
* [HPV genotyping](#hpv-genotyping) - reads alignment on all HPV strains
* [HPV local mapping](#hpv-local-mapping) - Local alignment on detected HPV strain(s)
* [Breakpoints detection](#breakpoints-detection) - Detection of soft-clipped ends on the HPV alignment file
* [Human mapping](#human-mapping) - Alignment of the soft-clipped reads to the Human reference genome
* [MultiQC](#multiqc) - aggregate report, describing results of the whole pipeline


## FastQC

Sequencing reads were first trimmed with [TrimGalore!](https://github.com/FelixKrueger/TrimGalore) to remove any adaptater sequences in the 3' end of the reads. This step is crutial to avoid noise in the downstream analysis and in the detection of soft-clipped reads.

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is then used to assess the overall sequencing quality. It provides information about the quality score distribution across your reads, the per base sequence content (%T/A/G/C). You get information about adapter contamination and other overrepresented sequences.

For further reading and documentation see the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

> **NB:** The FastQC plots displayed in the MultiQC report shows the trimmed reads.

**Output directory: `results/fastqc`**

* `sample_fastqc.html`
  * FastQC report, containing quality metrics for your untrimmed raw fastq files
* `zips/sample_fastqc.zip`
  * zip file containing the FastQC report, tab-delimited data file and plot images

## Control genes mapping

[Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) is used to map globally all the trimmed reads against 3 control genes (KLK3, GAPDH, RAB7A) in order to estimate virus load/burden.

For further reading and documentation see the [Bowtie2 manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml).

**Output directory: `results/ctrl_mapping`**

* `sample_ctrl.bam`
  * Mapping result of your trimmed fastq files against the 3 control genes
* `sample_ctrl.stats`
  * Statistics about the number of reads mapped per control gene

## HPV genotyping

[Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) is also used to map globally all the trimmed reads against all the HPVs genomes.

**Output directory: `results/hpv_mapping/global`**

* `sample_hpvs.bam`
  * Mapping result of your trimmed fastq files against all HPVs genomes
* `sample_HPVgenotyping.filtered`
  * List of at most 3 major hpv strains detected
* `sample_HPVgenotyping.stats`
  * Statistics about the number of reads mapped per HPV genome

## HPV local mapping

[Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) is finally used to map locally all the trimmed reads against the, at most 3, major HPVs strains previously detected genomes.

**Output directory: `results/hpv_mapping/local`**

* `sample-genotype.bam`
  * Mapping result of your trimmed fastq files against the corresponding genotype genome
* `sample-genotype_coverage.stats`
  * Table containing coverage statistic and composed of 5 columns : ID,sample,HPVsubtype,minCov,maxCov and meanCov
* `sample-genotype_covmatrix.mqc`
  * Coverage matrix (50bp bins)

## Breakpoints detection

**Output directory: `results/hpv_mapping/softclipped`**

* `sample-genotype_3prime_bkp.mqc`
  * This file contains for each position of the virus the corresponding number of reads having a 3' breakpoint at this position
* `sample-genotype_5prime_bkp.mqc`
  * This file contains for each position of the virus the corresponding number of reads having a 5' breakpoint at this position
* `sample-genotype.fa`
  * FASTA file containing all the human sequences of the hybrid reads with a detected breakpoint

## Human mapping

[Blat](https://genome.ucsc.edu/cgi-bin/hgBlat) is used to align all the human sequences of the hybrid reads with a detected breakpoint, on the human genome.

For further reading and documentation see the [Blat User Guide](https://genome.ucsc.edu/goldenPath/help/blatSpec.html).

**Output directory: `results/hpv_mapping/blat`**

* `sample-genotype.tsv`
  * Result of blat alignment
* `sample-genotype_table_to_display.csv`
  * Summary table that will be display in the MultiQC HTML report and containing breakpoint coordinates in both virus and human genome

## MultiQC

[MultiQC](http://multiqc.info) is a visualisation tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in within the report data directory.

The pipeline has special steps which allow the software versions used to be reported in the MultiQC output for future traceability.

**Output directory: `results/multiqc`**

* `Project_multiqc_report.html`
  * MultiQC report - a standalone HTML file that can be viewed in your web browser
* `Project_multiqc_data/`
  * Directory containing parsed statistics from the different tools used in the pipeline

For more information about how to use MultiQC reports, see http://multiqc.info
