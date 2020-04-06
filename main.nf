#!/usr/bin/env nextflow

/*
Copyright Institut Curie 2019
This software is a computer program whose purpose is to analyze high-throughput sequencing data.
You can use, modify and/ or redistribute the software under the terms of license (see the LICENSE file for more details).
The software is distributed in the hope that it will be useful, but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND. 
Users are therefore encouraged to test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data. 
The fact that you are presently reading this means that you have had knowledge of the license and that you accept its terms.

This script is based on the nf-core guidelines. See https://nf-co.re/ for more information
*/


/*
========================================================================================
                         HPV DETECTION PIPELINE
========================================================================================
 HPV Analysis Pipeline.
 #### Homepage / Documentation
 https://gitlab.curie.fr/data-analysis/illumina-hpv
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info"""
    
    HPV v${workflow.manifest.version}
    =======================================================

    Usage:

    nextflow run main.nf --reads '*_R{1,2}.fastq.gz' --genome 'hg19' -profile conda
    nextflow run main.nf --samplePlan sample_plan --genome hg19 -profile conda

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      --samplePlan                  Path to sample plan file if '--reads' is not specified
      --genome                      Name of iGenomes reference
      -profile                      Configuration profile to use. test / conda / toolsPath / singularity / cluster (see below)

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
      --nb_geno                     Number of HPV genotype to consider
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

     Skip options:
      --split_report                Generate one report per sample
      --skip_trimming               Skip trimming step
      --skip_fastqc                 Skip quality controls on sequencing reads
      --skip_blat                   Skip Human mapping with Blat
      --skip_multiqc                Skip report

    =======================================================                                                                                                                                                 
    Available Profiles
      -profile test                Set up the test dataset
      -profile conda               Build a new conda environment before running the pipeline
      -profile toolsPath           Use the paths defined in configuration for each tool
      -profile singularity         Use the Singularity images for each process
      -profile cluster             Run the workflow on the cluster, instead of locally

    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

// Configure reference genomes
// Reference index path configuration

params.bwt2_index = params.genome ? params.genomes[ params.genome ].bowtie2 ?: false : false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.blatdb =  params.genome ? params.genomes[ params.genome ].blatdb ?: false : false

params.bwt2_index_hpv = params.fasta_hpv ? false : params.genomes['HPV'].bowtie2 ?: false
params.bwt2_index_hpv_split = params.fasta_hpv ? false : params.genomes['HPV'].bowtie2_split ?: false
params.fasta_hpv = params.fasta_hpv ?: params.genomes['HPV'].fasta ?: false

params.genes_hpv = params.genomes['HPV'].genes ?: false
params.fasta_ctrl = params.genomes['HPV'].ctrl_capture ?: false



// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

// Stage config files
ch_multiqc_config = Channel.fromPath(params.multiqc_config)
ch_output_docs = Channel.fromPath("$baseDir/docs/output.md")
ch_genotype_nb = Channel.from(params.nb_geno)
ch_fasta_ctrl = Channel.fromPath(params.fasta_ctrl)
ch_hpv_genes_coord = Channel.fromPath(params.genes_hpv)

/*
 * CHANNELS
 */

/*
 * Create a channel for input read files
 */

if(params.samplePlan){
   if(params.singleEnd){
      Channel
         .from(file("${params.samplePlan}"))
         .splitCsv(header: false)
         .map{ row -> [ row[0], [file(row[2])]] }
         .set {reads_trimgalore}
   }else{
      Channel
         .from(file("${params.samplePlan}"))
         .splitCsv(header: false)
         .map{ row -> [ row[0], [file(row[2]), file(row[3])]] }
         .set {reads_trimgalore}
   }
   params.reads=false
}
else if(params.readPaths){
    if(params.singleEnd){
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .set {reads_trimgalore}
    } else {
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .set {reads_trimgalore}
    }
} else {
    Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
        .set {reads_trimgalore}
}

/*
 * Make sample plan if not available
 */

if (params.samplePlan){
  ch_splan = Channel.fromPath(params.samplePlan)
}else if(params.readPaths){
  if (params.singleEnd){
    Channel
       .from(params.readPaths)
       .collectFile() {
         item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + '\n']
        }
       .set{ ch_splan }
  }else{
     Channel
       .from(params.readPaths)
       .collectFile() {
         item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + ',' + item[1][1] + '\n']
        }
       .set{ ch_splan }
  }
}else{
  if (params.singleEnd){
    Channel
       .fromFilePairs( params.reads, size: 1 )
       .collectFile() {
          item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + '\n']
       }     
       .set { ch_splan }
  }else{
    Channel
       .fromFilePairs( params.reads, size: 2 )
       .collectFile() {
          item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + ',' + item[1][1] + '\n']
       }     
       .set { ch_splan }
   }
}

/*
 * Other input channels
 */

// Reference genome

if ( params.bwt2_index ){
   lastPath = params.bwt2_index.lastIndexOf(File.separator)
   ref_bwt2_dir =  params.bwt2_index.substring(0,lastPath+1)
   ref_bwt2_base = params.bwt2_index.substring(lastPath+1)

   Channel.fromPath( ref_bwt2_dir , checkIfExists: true)
      .ifEmpty { exit 1, "Genome index: Provided index not found: ${params.bwt2_index}" }
      .set { bwt2_ref_index }

}
else if ( params.fasta ) {
   lastPath = params.fasta.lastIndexOf(File.separator)
   ref_bwt2_base = params.fasta.substring(lastPath+1)

   Channel.fromPath( params.fasta )
        .ifEmpty { exit 1, "Genome index: Fasta file not found: ${params.fasta}" }
        .set { reference_fasta_for_index }
}
else {
   exit 1, "No reference genome specified!"
}

if ( params.blatdb ){
   Channel.fromPath( params.blatdb )
        .ifEmpty { exit 1, "BLAT database not found: ${params.blatdb}" }
        .set { blat_database }
}

//HPV genome

if ( params.bwt2_index_hpv && params.bwt2_index_hpv_split ){

   lastPath = params.bwt2_index_hpv.lastIndexOf(File.separator)
   hpv_bwt2_dir =  params.bwt2_index_hpv.substring(0,lastPath+1)
   hpv_bwt2_base = params.bwt2_index_hpv.substring(lastPath+1)

   Channel.fromPath( hpv_bwt2_dir , checkIfExists: true)
      .ifEmpty { exit 1, "HPV index: Provided index not found: ${params.bwt2_index_hpv}" }
      .set { bwt2_index_hpv }

   lastPath = params.bwt2_index_hpv_split.lastIndexOf(File.separator)
   hpv_split_bwt2_dir =  params.bwt2_index_hpv_split.substring(0,lastPath+1)
 
   Channel.fromPath( hpv_split_bwt2_dir , checkIfExists: true)
      .ifEmpty { exit 1, "HPV index per strain: Provided index not found: ${params.bwt2_index_hpv_split}" }
      .set { bwt2_index_hpv_split }
}
else if ( params.fasta_hpv ){
   lastPath = params.fasta_hpv.lastIndexOf(File.separator)
   hpv_bwt2_base = params.fasta_hpv.substring(lastPath+1) - ~/(\.fa)?(\.fasta)?(\.fas)?$/

   Channel.fromPath( params.fasta_hpv )
        .ifEmpty { exit 1, "HPV index: Fasta file not found: ${params.fasta_hpv}" }
        .set { hpv_fasta_for_index }
}
else{
   exit 1, "No HPV genome specified!"
}

// Header log info
log.info """=======================================================

HPV v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']  = 'HPV'
summary['Pipeline Version'] = workflow.manifest.version
summary['Run Name']     = custom_runName ?: workflow.runName
if (params.samplePlan) {
   summary['SamplePlan']   = params.samplePlan
}else{
   summary['Reads']        = params.reads
}
summary['Fasta Ref']    = params.fasta
summary['BLAT database']= params.blatdb
summary['Fasta HPV']    = params.fasta_hpv
summary['Split report'] = params.split_report
summary['Max Memory']   = params.max_memory
summary['Max CPUs']     = params.max_cpus
summary['Max Time']     = params.max_time
summary['Output dir']   = params.outdir
summary['Working dir']  = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
summary['Current user']   = "$USER"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Config Profile'] = workflow.profile

if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


/****************************************************
 * PRE-PROCESSING
 */

if ( !params.bwt2_index && params.fasta ){

  process makeBowtie2Index {
     tag "$ref_bwt2_base"
     publishDir path: { params.saveReference ? "${params.outdir}/references" : params.outdir },
              saveAs: { params.saveReference ? it : null }, mode: 'copy'
       
     input:
     file fasta from reference_fasta_for_index

     output:
     file "bowtie2_index" into bwt2_index_end2end
     file "bowtie2_index" into bwt2_index_trim

     script:
     ref_bwt2_base = fasta.toString() - ~/(\.fa)?(\.fasta)?(\.fas)?$/
     """
     mkdir bowtie2_index
     bowtie2-build ${fasta} bowtie2_index/${ref_bwt2_base}
     """
  }   
}

if ( (!params.bwt2_index_hpv | !params.bwt2_index_hpv_split) && params.fasta_hpv ){

  process makeBowtie2IndexHPV {
     publishDir path: { params.saveReference ? "${params.outdir}/references" : params.outdir },
              saveAs: { params.saveReference ? it : null }, mode: 'copy'
   
     input:
     file fasta from hpv_fasta_for_index

     output:
     file "bowtie2_index_hpv" into bwt2_index_hpv
     file "bowtie2_index_hpv_split" into bwt2_index_hpv_split

     script:
     hpv_bwt2_base = fasta.toString() - ~/(\.fa)?(\.fasta)?(\.fas)?$/
     """
     mkdir bowtie2_index_hpv
     bowtie2-build ${fasta} bowtie2_index_hpv/${hpv_bwt2_base}
   
     mkdir fasta_split && cd fasta_split
     cat ../$fasta | awk '{ if (substr(\$0, 1, 1)==">") {filename=(substr(\$0,2) ".fa")} print \$0 > filename }'
     cd .. && ls fasta_split/* > listoffasta.txt
   
     mkdir bowtie2_index_hpv_split
     while read ff; do
       base="\$(basename \"\${ff}\" | sed -e 's/.fa//')"
       bowtie2-build "\${ff}" bowtie2_index_hpv_split/"\${base}"
     done < listoffasta.txt
     """
  }
}


process makeBowtie2IndexCtrl {
  publishDir path: { params.saveReference ? "${params.outdir}/references" : params.outdir },
            saveAs: { params.saveReference ? it : null }, mode: 'copy'

  input:
  file fasta from ch_fasta_ctrl

  output:
  file "bowtie2_index_ctrl" into bwt2_index_ctrl

  script:
  """
  mkdir bowtie2_index_ctrl
  bowtie2-build ${fasta} bowtie2_index_ctrl/ctrl_regions
  """
}



/****************************************************
 * Main worflow
 */


/*
 * Reads Trimming
 */
if (!params.skip_trimming){
   process trimGalore {
     tag "$name" 

     publishDir "${params.outdir}/trimming", mode: 'copy',
                 saveAs: {filename -> filename.indexOf(".log") > 0 ? "logs/$filename" : "$filename"}

     input:
     set val(name), file(reads) from reads_trimgalore

     output:
     set val(name), file("*fq.gz") into reads_hpvmap, reads_splitmap, reads_ctrl, reads_fastqc
     set val(prefix), file("*trimming_report.txt") into trimgalore_results

     script:
     prefix = reads[0].toString() - ~/(_1)?(_2)?(_R1)?(_R2)?(.R1)?(.R2)?(_val_1)?(_val_2)?(\.fq)?(\.fastq)?(\.gz)?$/
     if (params.singleEnd) {
     """
     trim_galore --trim-n --quality 20 --length 20\
                --gzip $reads --basename ${prefix} --cores ${task.cpus}
     """
     }else {
     """
     trim_galore --trim-n --quality 20  --length 20 \
                --paired --gzip $reads --basename ${prefix} --cores ${task.cpus}
     mv ${prefix}_R1_val_1.fq.gz ${prefix}_R1_trimmed.fq.gz
     mv ${prefix}_R2_val_2.fq.gz ${prefix}_R2_trimmed.fq.gz
     mv ${reads[0]}_trimming_report.txt ${prefix}_R1_trimming_report.txt
     mv ${reads[1]}_trimming_report.txt ${prefix}_R2_trimming_report.txt
     """
     }
   }
}else{
   reads_trimgalore.into{reads_hpvmap; reads_splitmap; reads_ctrl; reads_fastqc}
   trimgalore_results = Channel.from(false)
}

/*
 * FastQC
 */
process fastqc {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}
   
    when:
    !params.skip_fastqc

    input:
    set val(name), file(reads) from reads_fastqc

    output:
    set val(prefix), file("${prefix}*.{zip,html}") into fastqc_results

    script:
    prefix = reads[0].toString() - ~/(_1)?(_2)?(_R1)?(_R2)?(.R1)?(.R2)?(_val_1)?(_val_2)?(_trimmed)?(\.fq)?(\.fastq)?(\.gz)?$/
    """
    fastqc -q $reads
    """
}

/*
 * Mapping
 */

process HPVmapping {
  tag "$prefix"
  publishDir "${params.outdir}/hpv_mapping/global", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".bam") == -1) "logs/$filename"
            else filename}

  input:
  set val(prefix), file(reads) from reads_hpvmap
  file index from bwt2_index_hpv.collect()

  output:
  set val(prefix), file("${prefix}_hpvs.bam") into hpv_bam

  script:
  if ( params.singleEnd ){
  """
  bowtie2 --rg-id BMG --rg SM:${prefix} \\
          --very-sensitive --no-unal \\
          -p ${task.cpus} \\
          -x ${index}/${hpv_bwt2_base} \\
          -U ${reads} > ${prefix}_hpvs.bam
  samtools sort -@  ${task.cpus} -o ${prefix}_sorted.bam ${prefix}_hpvs.bam
  mv ${prefix}_sorted.bam ${prefix}_hpvs.bam
  """
  }else{
  """
  bowtie2 --rg-id BMG --rg SM:${prefix} \\
          --very-sensitive --no-unal \\
          -p ${task.cpus} \\
          -x ${index}/${hpv_bwt2_base} \\
          -1 ${reads[0]} -2 ${reads[1]} > ${prefix}_hpvs.bam
  samtools sort -@  ${task.cpus} -o ${prefix}_sorted.bam ${prefix}_hpvs.bam
  mv ${prefix}_sorted.bam ${prefix}_hpvs.bam
  """
  }
}


process ctrlMapping {
  tag "$prefix"
  publishDir "${params.outdir}/ctrl_mapping/", mode: 'copy',
      saveAs: {filename ->
          if (filename.indexOf(".bam") == -1) "logs/$filename"
          else filename}

  input:
  set val(prefix), file(reads) from reads_ctrl
  file index from  bwt2_index_ctrl.collect()

  output:
  set val(prefix), file("${prefix}_ctrl.bam") into ctrl_bam

  script:
  if ( params.singleEnd ){
  """
  bowtie2 --rg-id BMG --rg SM:${prefix} \\
          --very-sensitive \\
          -p ${task.cpus} \\
          -x ${index}/ctrl_regions \\
          -U ${reads} > ${prefix}_ctrl.bam
  samtools sort -@  ${task.cpus} -o ${prefix}_sorted.bam ${prefix}_ctrl.bam
  mv ${prefix}_sorted.bam ${prefix}_ctrl.bam
  """
  }else{
  """
  bowtie2 --rg-id BMG --rg SM:${prefix} \\
          --very-sensitive \\
          -p ${task.cpus} \\
          -x ${index}/ctrl_regions \\
          -1 ${reads[0]} -2 ${reads[1]} > ${prefix}_ctrl.bam
  samtools sort -@  ${task.cpus} -o ${prefix}_sorted.bam ${prefix}_ctrl.bam
  mv ${prefix}_sorted.bam ${prefix}_ctrl.bam
  """
  }
}

process ctrlStats {
  tag "$prefix"
  publishDir "${params.outdir}/ctrl_mapping/", mode: 'copy'

  input:
  set val(prefix), file(bam) from ctrl_bam

  output:
  set val(prefix), file("*ctrl.stats") into ctrl_stats

  script:
  """
  samtools index ${bam}
  samtools idxstats ${bam} | cut -f1,3 | sort -k2,2nr > ${prefix}_ctrl.stats
  """
}

/*
 * Select Genotypes
 */

process selectGenotypes{
  publishDir "${params.outdir}/hpv_mapping/global", mode: 'copy'

  input:
  set val(prefix), file(bam) from hpv_bam 

  output:
  set val(prefix), file("${prefix}_HPVgenotyping.stats") into hpv_geno_stats
  file("${prefix}_HPVgenotyping.filtered") into sel_hpv_geno

  script:
  """
  samtools index ${bam}
  samtools idxstats ${bam} | cut -f1,3 | sort -k2,2nr > ${prefix}_HPVgenotyping.counts
  awk '\$2>=10{print}\$2<10{others+=\$2}END{print "Others\t"others}' ${prefix}_HPVgenotyping.counts | grep -v "*" > ${prefix}_HPVgenotyping.stats
  awk -F"[_\t]" '!(\$1 in done || \$1=="Others") {print \$0; done[\$1]=1}' ${prefix}_HPVgenotyping.stats | cut -f1 > ${prefix}_HPVgenotyping.filtered.tmp
  sed -i -e 's/^/${prefix},/' ${prefix}_HPVgenotyping.filtered.tmp
  head -n 3 ${prefix}_HPVgenotyping.filtered.tmp > ${prefix}_HPVgenotyping.filtered
  """
}


// Filter - removes all samples for which the genotype has not been detected
skipped_nogeno = []
def check_genotypes(geno) {
  def nb_geno = 0;
  geno.eachLine { nb_geno++; }
  samplename = geno.getBaseName() - '_HPVgenotyping.filered'
  if(nb_geno < 1 ){
      log.info "#################### NO HPV GENOTYPE DETECTED! IGNORING FOR FURTHER DOWNSTREAM ANALYSIS! ($samplename)"
      skipped_nogeno << samplename
      return false
  } else {
      return true
  }
}


sel_hpv_geno
        .filter { geno -> check_genotypes(geno) }
        .dump(tag: "toto")
    	.into { hpv_geno_filter; hpv_geno_mqc_config }

/*
 * Local Mapping for selected genotypes
 */

process HPVlocalMapping {

  publishDir "${params.outdir}/hpv_mapping/local", mode: 'copy',
      saveAs: {filename ->
          if (filename.indexOf(".bam") == -1) "logs/$filename"
          else filename}

  input:
  file index from bwt2_index_hpv_split.first()
  set val(prefix), val(hpv), file(reads) from hpv_geno_filter
    .splitCsv(header: ["sample", "hpv"])
    .map{
      [ it["sample"], it["hpv"] ]
    }
    .dump(tag: "mapped")
    .combine(reads_splitmap, by: 0)
    .dump(tag: "combined")

  output:
  set val(prefix), file("*.bam") into hpv_local_bam, hpv_cov_bam, hpv_soft_bam

  script: 
  if ( params.singleEnd ){
  """
  bowtie2 --rg-id BMG --rg SM:${prefix} \\
          --local --very-sensitive-local --no-unal \\
          -p ${task.cpus} \\
          -x ${index}/${hpv} \\
          -U ${reads} > ${prefix}-${hpv}.bam 2> ${prefix}-${hpv}_bowtie2.log
  samtools sort -@  ${task.cpus} -o ${prefix}-${hpv}_sorted.bam ${prefix}-${hpv}.bam
  mv ${prefix}-${hpv}_sorted.bam ${prefix}-${hpv}.bam
  """
  }else{
  """
  bowtie2 --rg-id BMG --rg SM:${prefix} \\
          --local --very-sensitive-local --no-unal \\
          -p ${task.cpus} \\
          -x ${index}/${hpv} \\
          -1 ${reads[0]} -2 ${reads[1]} > ${prefix}-${hpv}.bam 2> ${prefix}-${hpv}_bowtie2.log
  samtools sort -@  ${task.cpus} -o ${prefix}-${hpv}_sorted.bam ${prefix}-${hpv}.bam
  mv ${prefix}-${hpv}_sorted.bam ${prefix}-${hpv}.bam
  """
  }
}

process HPVlocalMappingStats {
  publishDir "${params.outdir}/hpv_mapping/local", mode: 'copy'
  
  input:
  set val(prefix), file(bam) from hpv_local_bam

  output:
  set val(prefix), file("*_coverage.stats") into hpv_cov_stats

  script:
  pfix= bam.toString() - ~/(.bam)?$/
  """
  genomeCoverageBed -d -ibam ${bam} > ${pfix}_coverage.out
  nbaln=\$(samtools flagstat ${bam} | grep 'mapped (' | awk '{print \$1}')
  echo -e "ID,sample,HPVsubtype,mappedReads,minCov,maxCov,meanCov" > ${pfix}_coverage.stats
  awk -v id=${prefix} -v nbaln=\$nbaln 'NR==1{hpv=\$1;mn=mx=\$3}{total+=\$3}(\$3>mx){mx=\$3}(\$3<mn){mn=\$3} END{OFS=","; print id"_"hpv,id,hpv,nbaln,mn,mx,total/NR}' ${pfix}_coverage.out >> ${pfix}_coverage.stats
  """
}


process HPVcoverage {
  publishDir "${params.outdir}/hpv_mapping/local", mode: 'copy'

  input:
  set val(prefix), file(bam) from hpv_cov_bam

  output:
  set val(prefix), file("*covmatrix.mqc") into hpv_bw_cov

  script:
  pfix= bam.toString() - ~/(_sorted)?(.bam)?$/
  """
  samtools index ${bam}
  bamCoverage -b ${bam} --binSize 50 --normalizeUsing CPM --outFileFormat bedgraph -o ${pfix}.bedgraph
  awk -F"\t" '{OFS="\t"; print \$2+25,\$4}' ${pfix}.bedgraph > ${pfix}_covmatrix.mqc
  """
}


/*
 * Breakpoint detection
 */

process extractBreakpointsSequence {
   publishDir "${params.outdir}/hpv_mapping/softclipped", mode: 'copy'

   input:
   set val(prefix), file(bam) from hpv_soft_bam

   output:
   set val(prefix), file("*.mqc") into bkp_pos mode 'flatten'
   set val(prefix), file("*.fa") into clipped_seq

   script:
   pfix= bam.toString() - ~/.bam$/
   """
   extractSoftclipped.py -v --mqc --stranded --minLen 10 ${bam}
   sort -k1,1n ${pfix}_3prime_bkp.mqc | awk 'BEGIN{nr=1} nr<\$1{for (i=nr;i<\$1;i++){print i"\t"0} nr=\$1}{print; nr+=1}' > file1.tmp
   mv file1.tmp ${pfix}_3prime_bkp.mqc
   sort -k1,1n ${pfix}_5prime_bkp.mqc | awk 'BEGIN{nr=1} nr<\$1{for (i=nr;i<\$1;i++){print i"\t"0} nr=\$1}{print; nr+=1}' > file2.tmp
   mv file2.tmp ${pfix}_5prime_bkp.mqc
   """
}

/*
 * Blat
 */  

if (!params.skip_blat){
   process blatSoftClippedSeq {
      publishDir "${params.outdir}/hpv_mapping/blat", mode: 'copy'

      when:
      !params.skip_blat

      input:
      file(blatdb) from blat_database.collect()
      set val(prefix),file(fasta) from clipped_seq
 
      output:
      set val(prefix), file("*.tsv") into blat_res

      script:
      pfix= fasta.toString() - ~/(.fa)?$/
      """
      blat ${blatdb} ${fasta} ${pfix}.tsv -noHead -minScore=25 -minIdentity=90
      """
   }

   process blatSummary {
      publishDir "${params.outdir}/hpv_mapping/blat", mode: 'copy'

      input:
      set val(prefix), file(psl) from blat_res
 
      output:
      set val(prefix), file("*_table_filtered.csv") into ttd
      set val(prefix), file("*_table.csv") into table

      script:
      pfix= psl.toString() - ~/(.tsv)?$/
      """
      blat_parser.py ${psl}
      """
   }
}else{
   ttd = Channel.from(false)
}


/*
/* MultiQC
 */


process get_software_versions {

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    bowtie2 --version > v_bowtie2.txt
    samtools --version > v_samtools.txt
    bedtools --version > v_bedtools.txt
    deeptools --version 2> v_deeptools.txt
    echo "BLAT v. 35" > v_blat.txt
    python --version 2> v_python.txt
    multiqc --version > v_multiqc.txt
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}

process workflow_summary_mqc {
  when:
  !params.skip_multiqc

  output:
  file 'workflow_summary_mqc.yaml' into workflow_summary_yaml

  exec:
  def yaml_file = task.workDir.resolve('workflow_summary_mqc.yaml')
  yaml_file.text  = """
  id: 'summary'
  description: " - this information is collected when the pipeline is started."
  section_name: 'Workflow Summary'
  section_href: 'https://gitlab.curie.fr/illumina-hpv'
  plot_type: 'html'
  data: |
      <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
      </dl>
  """.stripIndent()
}

if (params.split_report){

   process make_hpv_config_persample {        
      when:
      !params.skip_multiqc

      input:
      file(geno) from hpv_geno_mqc_config
      file genes from ch_hpv_genes_coord.collect()

      output:
      set val(prefix),file('*conf.mqc') into mqc_hpv_conf
      set val(prefix),file('*bkp.mqc') into mqc_genepos

      script:
      prefix = geno.toString() - "_HPVgenotyping.filtered"
      """
      awk -F, '{print \$2}' ${geno} | sort -u > allgenotypes_unique.txt
      scrape_mqc_config.py  allgenotypes_unique.txt > ${prefix}_hpv_conf.mqc
      gene_tracks.sh allgenotypes_unique.txt ${genes} ${prefix} 
      """
   }

   mqc_hpv_conf
        .join(fastqc_results, remainder: true)
        .join(trimgalore_results, remainder: true)
        .join(ctrl_stats)
        .join(hpv_cov_stats.groupTuple(), remainder: true)
        .join(hpv_bw_cov.groupTuple(), remainder: true)
        .join(bkp_pos.groupTuple(), remainder: true)
        .join(hpv_geno_stats, remainder: true)
        .join(mqc_genepos, remainder: true)
	.join(ttd.groupTuple(), remainder: true)
	.dump(tag: "join")
        .set{ch_hpv_report}

   process multiqc {
     publishDir "${params.outdir}/MultiQC/", mode: 'copy'

     when:
     !params.skip_multiqc

     input:
     file splan from ch_splan.first()
     file multiqc_config from ch_multiqc_config.first()
     set val(prefix), file('mqc/hpv_config.mqc'), file('fastqc/*'), file('trimming/*'), file('ctrl/*'), file('hpv/*'), file('hpv/*'), file('hpv/*'), file('hpv/*'), file('mqc/*'), file('mqc/*') from ch_hpv_report.dump(tag: "mqc") 
     file ('software_versions/*') from software_versions_yaml.collect()
     file ('workflow_summary/*') from workflow_summary_yaml.collect()

     output:
     file splan
     file "*.html" into multiqc_report
     file "*_data"

     script:
     rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
     rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_" + prefix + "_multiqc_report" : ''
     """	
     multiqc . -f $rtitle -n ${prefix}_HPVreport.html -c $multiqc_config -c 'mqc/hpv_config.mqc' -m fastqc -m custom_content
     """
   }
}else{
   process make_hpv_config {

      when:
      !params.skip_multiqc

      input:
      file(geno) from hpv_geno_mqc_config.collect()
      file genes from ch_hpv_genes_coord.collect()

      output:
      file('*conf.mqc') into mqc_hpv_conf
      file('*bkp.mqc') into mqc_genepos

      script:
      prefix = geno.toString() - "_HPVgenotyping.filtered"
      """
      awk -F, '{print \$2}' ${geno} | sort -u > allgenotypes_unique.txt
      scrape_mqc_config.py  allgenotypes_unique.txt > hpv_conf.mqc
      gene_tracks.sh allgenotypes_unique.txt ${genes} 'genes' 
      """
   }

   process multiqc_allsamples {
     publishDir "${params.outdir}/MultiQC/", mode: 'copy'

     when:
     !params.skip_multiqc

     input:
     file splan from ch_splan.first()
     file multiqc_config from ch_multiqc_config.first()
     file hpv_config from mqc_hpv_conf.collect().ifEmpty([])
     file('fastqc/*') from fastqc_results.map{items->items[1]}.collect().ifEmpty([])
     file('trimming/*') from trimgalore_results.map{items->items[1]}.collect().ifEmpty([])
     file ('hpv/*') from hpv_geno_stats.map{items->items[1]}.collect()
     file ('hpv/*') from hpv_cov_stats.map{items->items[1]}.collect()
     file ('hpv/*') from hpv_bw_cov.map{items->items[1]}.collect()
     file ('hpv/*') from bkp_pos.map{items->items[1]}.collect()
     file ('hpv/*') from mqc_genepos.collect()
     file ('hpv/*') from ttd.collect().ifEmpty([])
     file ('ctrl/*') from ctrl_stats.collect()
  
     file ('software_versions/*') from software_versions_yaml.collect()
     file ('workflow_summary/*') from workflow_summary_yaml.collect()
 
     output:
     file splan
     file "*multiqc_report.html" into multiqc_report
     file "*_data"

     script:
     rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
     rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
     metadata_opts = params.metadata ? "--metadata ${metadata}" : ""
     splan_opts = params.samplePlan ? "--splan ${params.samplePlan}" : ""
     """	
     mqc_header.py --name "nv-VIF" --version ${workflow.manifest.version} ${metadata_opts} ${splan_opts} > multiqc-config-header.yaml         
     multiqc . -f $rtitle $rfilename -c $multiqc_config -c $hpv_config -c multiqc-config-header.yaml -m fastqc -m custom_content
     """
   }
}
