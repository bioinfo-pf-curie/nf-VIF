#!/usr/bin/env nextflow

/*
Copyright Institut Curie 2019-2020
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
      --bwt2Index                   Path to Bowtie2 index
      --fasta                       Path to Fasta reference (.fasta)
      --blatdb                      Path to BLAT database (.2bit)

    HPV References:
      --fastaHpv                    Path to Fasta HPV reference (.fasta)                 
      --bwt2IndexHpv                Path to Bowtie2 index for all HPV strains
      --bwt2IndexHpvSplit           Path to Bowtie2 index per HPV strain
      --saveReference               Save all references generated during the analysis. Default: False

    Advanced options:
      --minMapq                     Minimum reads mapping quality. Default: 0
      --minLen                      Minimum trimmed length sequence to consider. Default: 15
      --minFreqGeno                 Fraction of reads to consider a genotpye. Default: 0.2
      --nbGeno                      Number of HPV genotype to consider
      --splitReport                 Generate one report per sample
 
    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

     Skip options:
      --skipTrimming               Skip trimming step
      --skipFastqc                 Skip quality controls on sequencing reads
      --skipBlat                   Skip Human mapping with Blat
      --skipMultiqc                Skip report

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

params.bwt2Index = params.genome ? params.genomes[ params.genome ].bowtie2 ?: false : false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.blatdb =  params.genome ? params.genomes[ params.genome ].blatdb ?: false : false

params.bwt2IndexHpv = params.fastaHpv ? false : params.genomes['HPV'].bowtie2 ?: false
params.bwt2IndexHpvSplit = params.fastaHpv ? false : params.genomes['HPV'].bowtie2Split ?: false
params.fastaHpv = params.fastaHpv ?: params.genomes['HPV'].fasta ?: false

params.genesHpv = params.genomes['HPV'].genes ?: false
params.fastaCtrl = params.genomes['HPV'].ctrlCapture ?: false



// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
customRunName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  customRunName = workflow.runName
}

// Stage config files
chMultiqcConfig = Channel.fromPath(params.multiqcConfig)
chOutputDocs = Channel.fromPath("$baseDir/docs/output.md")
chFastaCtrl = Channel.fromPath(params.fastaCtrl)
chHpvGenesCoord = Channel.fromPath(params.genesHpv)

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
         .set {readsTrimgalore}
   }else{
      Channel
         .from(file("${params.samplePlan}"))
         .splitCsv(header: false)
         .map{ row -> [ row[0], [file(row[2]), file(row[3])]] }
         .set {readsTrimgalore}
   }
   params.reads=false
}
else if(params.readPaths){
    if(params.singleEnd){
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .set {readsTrimgalore}
    } else {
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .set {readsTrimgalore}
    }
} else {
    Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
        .set {readsTrimgalore}
}

/*
 * Make sample plan if not available
 */

if (params.samplePlan){
  chSplan = Channel.fromPath(params.samplePlan)
}else if(params.readPaths){
  if (params.singleEnd){
    Channel
       .from(params.readPaths)
       .collectFile() {
         item -> ["samplePlan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + '\n']
        }
       .set{ chSplan }
  }else{
     Channel
       .from(params.readPaths)
       .collectFile() {
         item -> ["samplePlan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + ',' + item[1][1] + '\n']
        }
       .set{ chSplan }
  }
}else{
  if (params.singleEnd){
    Channel
       .fromFilePairs( params.reads, size: 1 )
       .collectFile() {
          item -> ["samplePlan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + '\n']
       }     
       .set { chSplan }
  }else{
    Channel
       .fromFilePairs( params.reads, size: 2 )
       .collectFile() {
          item -> ["samplePlan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + ',' + item[1][1] + '\n']
       }     
       .set { chSplan }
   }
}

/*
 * Other input channels
 */

// Reference genome

if ( params.bwt2Index ){
   lastPath = params.bwt2Index.lastIndexOf(File.separator)
   refBwt2Dir =  params.bwt2Index.substring(0,lastPath+1)
   refBwt2Base = params.bwt2Index.substring(lastPath+1)

   Channel.fromPath( refBwt2Dir , checkIfExists: true)
      .ifEmpty { exit 1, "Genome index: Provided index not found: ${params.bwt2Index}" }
      .set { bwt2RefIndex }

}
else if ( params.fasta ) {
   lastPath = params.fasta.lastIndexOf(File.separator)
   refBwt2Base = params.fasta.substring(lastPath+1)

   Channel.fromPath( params.fasta )
        .ifEmpty { exit 1, "Genome index: Fasta file not found: ${params.fasta}" }
        .set { referenceFastaForIndex }
}
else {
   exit 1, "No reference genome specified!"
}

if ( params.blatdb ){
   Channel.fromPath( params.blatdb )
        .ifEmpty { exit 1, "BLAT database not found: ${params.blatdb}" }
        .set { blatDatabase }
}

//HPV genome

if ( params.bwt2IndexHpv && params.bwt2IndexHpvSplit ){

   lastPath = params.bwt2IndexHpv.lastIndexOf(File.separator)
   hpvBwt2Dir =  params.bwt2IndexHpv.substring(0,lastPath+1)
   hpvBwt2Base = params.bwt2IndexHpv.substring(lastPath+1)

   Channel.fromPath( hpvBwt2Dir , checkIfExists: true)
      .ifEmpty { exit 1, "HPV index: Provided index not found: ${params.bwt2IndexHpv}" }
      .set { bwt2IndexHpv }

   lastPath = params.bwt2IndexHpvSplit.lastIndexOf(File.separator)
   hpvSplitBwt2Dir =  params.bwt2IndexHpvSplit.substring(0,lastPath+1)
 
   Channel.fromPath( hpvSplitBwt2Dir , checkIfExists: true)
      .ifEmpty { exit 1, "HPV index per strain: Provided index not found: ${params.bwt2IndexHpvSplit}" }
      .set { bwt2IndexHpvSplit }
}
else if ( params.fastaHpv ){
   lastPath = params.fastaHpv.lastIndexOf(File.separator)
   hpvBwt2Base = params.fastaHpv.substring(lastPath+1) - ~/(\.fa)?(\.fasta)?(\.fas)?$/

   Channel.fromPath( params.fastaHpv )
        .ifEmpty { exit 1, "HPV index: Fasta file not found: ${params.fastaHpv}" }
        .set { hpvFastaForIndex }
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
summary['Run Name']     = customRunName ?: workflow.runName
if (params.samplePlan) {
   summary['SamplePlan']   = params.samplePlan
}else{
   summary['Reads']        = params.reads
}
summary['Fasta Ref']      = params.fasta
summary['BLAT database']  = params.blatdb
summary['Fasta HPV']      = params.fastaHpv
summary['Min MAPQ']       = params.minMapq
summary['Min length']     = params.minLen
summary['Min Freq Geno']  = params.minFreqGeno
summary['Split report']   = params.splitReport
summary['Max Memory']     = params.max_memory
summary['Max CPUs']       = params.max_cpus
summary['Max Time']       = params.max_time
summary['Output dir']     = params.outdir
summary['Working dir']    = workflow.workDir
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

if ( !params.bwt2Index && params.fasta ){

  process makeBowtie2Index {
     tag "$refBwt2Base"
     publishDir path: { params.saveReference ? "${params.outdir}/references" : params.outdir },
              saveAs: { params.saveReference ? it : null }, mode: 'copy'
       
     input:
     file fasta from referenceFastaForIndex

     output:
     file "bowtie2Index" into bwt2IndexEnd2end
     file "bowtie2Index" into bwt2IndexTrim

     script:
     refBwt2Base = fasta.toString() - ~/(\.fa)?(\.fasta)?(\.fas)?$/
     """
     mkdir bowtie2Index
     bowtie2-build ${fasta} bowtie2Index/${refBwt2Base}
     """
  }   
}

if ( (!params.bwt2IndexHpv | !params.bwt2IndexHpvSplit) && params.fastaHpv ){

process makeBowtie2IndexHPV {
     publishDir path: { params.saveReference ? "${params.outdir}/references" : params.outdir },
              saveAs: { params.saveReference ? it : null }, mode: 'copy'
   
     input:
     file fasta from hpvFastaForIndex

     output:
     file "bowtie2IndexHpv" into bwt2IndexHpv
     file "bowtie2IndexHpvSplit" into bwt2IndexHpvSplit

     script:
     hpvBwt2Base = fasta.toString() - ~/(\.fa)?(\.fasta)?(\.fas)?$/
     """
     mkdir bowtie2IndexHpv
     bowtie2-build ${fasta} bowtie2IndexHpv/${hpvBwt2Base}
   
     mkdir fastaSplit && cd fastaSplit
     cat ../$fasta | awk '{ if (substr(\$0, 1, 1)==">") {filename=(substr(\$0,2) ".fa")} print \$0 > filename }'
     cd .. && ls fastaSplit/* > listoffasta.txt
   
     mkdir bowtie2IndexHpvSplit
     while read ff; do
       base="\$(basename \"\${ff}\" | sed -e 's/.fa//')"
       bowtie2-build "\${ff}" bowtie2IndexHpvSplit/"\${base}"
     done < listoffasta.txt
     """
}
}


process makeBowtie2IndexCtrl {
  publishDir path: { params.saveReference ? "${params.outdir}/references" : params.outdir },
            saveAs: { params.saveReference ? it : null }, mode: 'copy'

  input:
  file fasta from chFastaCtrl

  output:
  file "bowtie2IndexCtrl" into bwt2IndexCtrl

  script:
  """
  mkdir bowtie2IndexCtrl
  bowtie2-build ${fasta} bowtie2IndexCtrl/ctrlRegions
  """
}



/****************************************************
 * Main worflow
 */


/*
 * Reads Trimming
 */
if (!params.skipTrimming){
   process trimGalore {
     tag "$name" 

     publishDir "${params.outdir}/trimming", mode: 'copy',
                 saveAs: {filename -> filename.indexOf(".log") > 0 ? "logs/$filename" : "$filename"}

     input:
     set val(name), file(reads) from readsTrimgalore

     output:
     set val(name), file("*fq.gz") into readsHpvmap, readsSplitmap, readsCtrl, readsFastqc
     set val(prefix), file("*trimming_report.txt") into trimgaloreResults

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
   readsTrimgalore.into{readsHpvmap; readsSplitmap; readsCtrl; readsFastqc}
   trimgaloreResults = Channel.from(false)
}


/*
 * FastQC
 */
process fastqc {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}
   
    when:
    !params.skipFastqc

    input:
    set val(name), file(reads) from readsFastqc

    output:
    set val(prefix), file("${prefix}*.{zip,html}") into fastqcResults

    script:
    prefix = reads[0].toString() - ~/(_1)?(_2)?(_R1)?(_R2)?(.R1)?(.R2)?(_val_1)?(_val_2)?(_trimmed)?(\.fq)?(\.fastq)?(\.gz)?$/
    """
    fastqc -q $reads
    """
}

/*
 * Mapping on control regions
 */

process ctrlMapping {
  tag "$prefix"
  publishDir "${params.outdir}/ctrlMapping/", mode: 'copy',
      saveAs: {filename ->
          if (filename.endsWith(".log")) "logs/$filename"
          else if (params.saveAlignedIntermediates) filename
	  else null}

  input:
  set val(prefix), file(reads) from readsCtrl
  file index from  bwt2IndexCtrl.collect()

  output:
  set val(prefix), file("${prefix}_ctrl.bam") into ctrlBam
  set val(prefix), file("*ctrl_bowtie2.log") into ctrlBowtie2Log

  script:
  if ( params.singleEnd ){
  """
  bowtie2 --rg-id BMG --rg SM:${prefix} \\
          --very-sensitive \\
          -p ${task.cpus} \\
          -x ${index}/ctrlRegions \\
          -U ${reads} > ${prefix}_ctrl.bam 2> ${prefix}_ctrl_bowtie2.log
  """
  }else{
  """
  bowtie2 --rg-id BMG --rg SM:${prefix} \\
          --very-sensitive \\
          -p ${task.cpus} \\
          -x ${index}/ctrlRegions \\
          -1 ${reads[0]} -2 ${reads[1]} > ${prefix}_ctrl.bam 2> ${prefix}_ctrl_bowtie2.log
  """
  }
}

process ctrlStats {
  tag "$prefix"
  publishDir "${params.outdir}/ctrlMapping/", mode: 'copy'

  input:
  set val(prefix), file(bam) from ctrlBam

  output:
  set val(prefix), file('*fsorted_ctrl.{bam,bam.bai}') into ctrlFiltBams
  set val(prefix), file("*ctrl.stats") into ctrlStats

  script:
  peStatus = params.singleEnd ? "0" : "1"
  """
  nbreads=\$(samtools view -c ${bam})
  samtools view -h -q 20 ${bam} | samtools sort -@  ${task.cpus} - > ${prefix}_fsorted_ctrl.bam
  samtools index ${prefix}_fsorted_ctrl.bam
  samtools idxstats ${prefix}_fsorted_ctrl.bam | cut -f1,3 | sort -k2,2nr > ${prefix}_ctrl.stats
  awk -v isPe=${peStatus} -v tot=\$nbreads '\$1!="*"{s=s+\$2} \$1=="*"{\$1="unmapped"; \$2=tot-s} isPe==1{\$2=\$2/2} {printf "%s\\t%.0f\\n",\$1,\$2}' ${prefix}_ctrl.stats > ${prefix}_ctrl_final.stats
  mv ${prefix}_ctrl_final.stats ${prefix}_ctrl.stats
  """
}


/*
 * HPV mapping and genotyping
 */ 

process HPVmapping {
  tag "$prefix"
  publishDir "${params.outdir}/hpvMapping/allref", mode: 'copy',
        saveAs: {filename ->
            if (filename.endsWith(".log")) "logs/$filename"
            else if (params.saveAlignedIntermediates) filename
	    else null}

  input:
  set val(prefix), file(reads) from readsHpvmap
  file index from bwt2IndexHpv.collect()

  output:
  set val(prefix), file("${prefix}_hpvs.bam") into hpvBam
  set val(prefix), file("*bowtie2.log") into hpvBowtie2Log

  script:
  if ( params.singleEnd ){
  """
  bowtie2 --rg-id BMG --rg SM:${prefix} \\
          --very-sensitive --no-unal \\
          -p ${task.cpus} \\
          -x ${index}/${hpvBwt2Base} \\
          -U ${reads} > ${prefix}_hpvs.bam 2> ${prefix}_hpvs_bowtie2.log
  """
  }else{
  """
  bowtie2 --rg-id BMG --rg SM:${prefix} \\
          --very-sensitive --no-unal \\
          -p ${task.cpus} \\
          -x ${index}/${hpvBwt2Base} \\
          -1 ${reads[0]} -2 ${reads[1]} > ${prefix}_hpvs.bam 2> ${prefix}_hpvs_bowtie2.log 
  """
  }
}

process selectGenotypes{
  publishDir "${params.outdir}/hpvMapping/allref", mode: 'copy'

  input:
  set val(prefix), file(bam) from hpvBam 

  output:
  set val(prefix), file("${prefix}_HPVgenotyping.stats") into hpvGenoStats
  set val(prefix), file("${prefix}_HPVgenotyping.filtered") into hpvGenoMqc
  file("${prefix}_HPVgenotyping.filtered") into selHpvGeno
  set val(prefix), file('*fsorted_hpvs.{bam,bam.bai}') into hpvsFiltBams

  script:
  """
  samtools view -h -q ${params.minMapq} ${bam} | samtools sort -@  ${task.cpus} -o ${prefix}_fsorted_hpvs.bam -
  samtools index ${prefix}_fsorted_hpvs.bam
  samtools idxstats ${prefix}_fsorted_hpvs.bam | cut -f1,3 | sort -k2,2nr > ${prefix}_HPVgenotyping.counts
  nbreads=\$(samtools view -c ${bam})
  awk -v tot=\$nbreads '{printf("%s\\t%.02f\\n", \$1, \$2/tot)}' ${prefix}_HPVgenotyping.counts > ${prefix}_HPVgenotyping.freq
  awk -v minFreq=${params.minFreqGeno} -v sname=${prefix} '\$2>=minFreq{print sname","\$1}' ${prefix}_HPVgenotyping.freq > ${prefix}_HPVgenotyping.filtered
  awk '\$2>=10{print}\$2<10{others+=\$2}END{print "Others\t"others}' ${prefix}_HPVgenotyping.counts | grep -v "*" > ${prefix}_HPVgenotyping.stats
  """
}


// Filter - removes all samples for which the genotype has not been detected
skippedNogeno = []
def checkGenotypes(geno) {
  def nbGeno = 0;
  geno.eachLine { nbGeno++; }
  samplename = geno.getBaseName() - '_HPVgenotyping.filered'
  if(nbGeno < 1 ){
      log.info "#################### NO HPV GENOTYPE DETECTED! IGNORING FOR FURTHER DOWNSTREAM ANALYSIS! ($samplename)"
      skippedNogeno << samplename
      return false
  } else {
      return true
  }
}


selHpvGeno
        .filter { geno -> checkGenotypes(geno) }
    	.into { hpvGenoFilter; hpvGenoMqcConfig }

/*
 * Local Mapping for selected genotypes
 */

process HPVlocalMapping {
  publishDir "${params.outdir}/hpvMapping/pergenotype", mode: 'copy',
      saveAs: {filename ->
          if (filename.endsWith(".log")) "logs/$filename"
	  else if (params.saveAlignedIntermediates) filename
          else null}

  input:
  file index from bwt2IndexHpvSplit.first()
  set val(prefix), val(hpv), file(reads) from hpvGenoFilter
    .splitCsv(header: ["sample", "hpv"])
    .map{
      [ it["sample"], it["hpv"] ]
    }
    .combine(readsSplitmap, by: 0)
    .dump(tag: "hpvloc")

  output:
  set val(prefix), file("*.bam") into hpvLocalBam, hpvCovBam, hpvSoftBam

  script: 
  if ( params.singleEnd ){
  """
  bowtie2 --rg-id BMG --rg SM:${prefix} \\
          --local --very-sensitive-local --no-unal \\
          -p ${task.cpus} \\
          -x ${index}/${hpv} \\
          -U ${reads} > ${prefix}-${hpv}.bam 2> ${prefix}-${hpv}_bowtie2.log
  """
  }else{
  """
  bowtie2 --rg-id BMG --rg SM:${prefix} \\
          --local --very-sensitive-local --no-unal \\
          -p ${task.cpus} \\
          -x ${index}/${hpv} \\
          -1 ${reads[0]} -2 ${reads[1]} > ${prefix}-${hpv}.bam 2> ${prefix}-${hpv}_bowtie2.log
  """
  }
}


process HPVlocalMappingStats {
  publishDir "${params.outdir}/hpvMapping/pergenotype", mode: 'copy'
  
  input:
  set val(prefix), file(bam) from hpvLocalBam

  output:
  set val(prefix), file("*_coverage.stats") into hpvCovStats

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
  publishDir "${params.outdir}/hpvMapping/pergenotype", mode: 'copy'

  input:
  set val(prefix), file(bam) from hpvCovBam

  output:
  set val(prefix), file("*covmatrix.mqc") into hpvBwCov
  set val(prefix), file('*sorted.{bam,bam.bai}') into hpvSortedBams

  script:
  pfix= bam.toString() - ~/(_sorted)?(.bam)?$/
  normOpts = params.splitReport ? "" : "--normalizeUsing CPM"
  """
  samtools sort -@ ${task.cpus} -o ${pfix}_sorted.bam ${bam}
  samtools index ${pfix}_sorted.bam
  bamCoverage -b ${pfix}_sorted.bam --binSize 50 ${normOpts} --outFileFormat bedgraph -o ${pfix}.bedgraph
  awk -F"\t" '{OFS="\t"; print \$2+25,\$4}' ${pfix}.bedgraph > ${pfix}_covmatrix.mqc
  """
}


/*
 * Breakpoint detection
 */

process extractBreakpointsSequence {
   publishDir "${params.outdir}/hpvMapping/softclipped", mode: 'copy',
               saveAs: {filename -> 
                   if (filename.indexOf(".mqc") > 0) "mqc/$filename"
		   else filename}

   input:
   set val(prefix), file(bam) from hpvSoftBam

   output:
   set val(prefix), file("*.mqc") into bkpPos mode 'flatten'
   set val(pfix), file("*.csv") into bkpInfo
   set val(pfix), val(prefix), file("*.fa") into clippedSeq

   script:
   pfix= bam.toString() - ~/.bam$/
   """
   extractSoftclipped.py -v --mqc --stranded --minLen ${params.minLen} ${bam}
   sort -k1,1n ${pfix}_3prime_bkp.mqc | awk 'BEGIN{nr=1} nr<\$1{for (i=nr;i<\$1;i++){print i"\t"0} nr=\$1}{print; nr+=1}' > file1.tmp
   mv file1.tmp ${pfix}_3prime_bkp.mqc
   sort -k1,1n ${pfix}_5prime_bkp.mqc | awk 'BEGIN{nr=1} nr<\$1{for (i=nr;i<\$1;i++){print i"\t"0} nr=\$1}{print; nr+=1}' > file2.tmp
   mv file2.tmp ${pfix}_5prime_bkp.mqc
   """
}

/*
 * Blat
 */  

if (!params.skipBlat){
   process blatSoftClippedSeq {
      publishDir "${params.outdir}/hpvMapping/blat", mode: 'copy'

      when:
      !params.skipBlat

      input:
      file(blatdb) from blatDatabase.collect()
      set val(pfix), val(sname), file(fasta) from clippedSeq
 
      output:
      set val(pfix), val(sname), file("*.tsv") into blatRes

      script:
      """
      blat ${blatdb} ${fasta} ${pfix}.tsv -noHead -minScore=25 -minIdentity=90
      """
   }
   
   process blatSummary {
      publishDir "${params.outdir}/hpvMapping/blat", mode: 'copy'

      input:
      set val(pfix), val(sname), file(psl), file(csv) from blatRes.join(bkpInfo).dump(tag:"blat")
 
      output:
      set val(sname), file("*_table_filtered.csv") into ttdHQ
      set val(sname), file("*_table.csv") into tableHQ
      set val(sname), file("*_bkptable_filtered.csv") into ttd
      set val(sname), file("*_bkptable.csv") into table

      script:
      """
      blatParser.py -f ${psl} -b ${csv} --sname ${sname}
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
  !params.skipMultiqc

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

if (params.splitReport){

   process makeHpvConfigPerSample {        
      when:
      !params.skipMultiqc

      input:
      file(geno) from hpvGenoMqcConfig
      file genes from chHpvGenesCoord.collect()

      output:
      set val(prefix),file('*conf.mqc') into mqcHpvConf
      set val(prefix),file('*bkp.mqc') optional true into mqcGenepos

      script:
      prefix = geno.toString() - "_HPVgenotyping.filtered"
      """
      awk -F, '{print \$2}' ${geno} | sort -u > allgenotypes_unique.txt
      scrape_mqc_config.py  allgenotypes_unique.txt > ${prefix}_hpv_conf.mqc
      gene_tracks.sh allgenotypes_unique.txt ${genes} ${prefix} 
      """
   }

   mqcHpvConf
        .join(fastqcResults, remainder: true)
        .join(trimgaloreResults, remainder: true)
        .join(ctrlStats)
	.join(hpvBowtie2Log)
	.join(hpvGenoMqc)
        .join(hpvCovStats.groupTuple(), remainder: true)
        .join(hpvBwCov.groupTuple(), remainder: true)
        .join(bkpPos.groupTuple(), remainder: true)
        .join(hpvGenoStats, remainder: true)
        .join(mqcGenepos, remainder: true)
	.join(ttd.groupTuple(), remainder: true)
	.dump(tag: "join")
        .set{chHpvReport}

   process multiqc {
     publishDir "${params.outdir}/MultiQC/", mode: 'copy'

     when:
     !params.skipMultiqc

     input:
     file splan from chSplan.first()
     file multiqcConfig from chMultiqcConfig.first()
     set val(prefix), file('mqc/hpv_config.mqc'), file('fastqc/*'), file('trimming/*'), file('ctrl/*'), 
     file('hpv/*'), file('hpv/*'), file('hpv/*'), file('hpv/*'), file('hpv/*'), file('hpv/*'), file('hpv/*'), file('hpv/*') from chHpvReport.dump(tag: "mqc") 
     file ('software_versions/*') from software_versions_yaml.collect()
     file ('workflow_summary/*') from workflow_summary_yaml.collect()

     output:
     file splan
     file "*.html" into multiqc_report
     file "*_data"

     script:
     rtitle = customRunName ? "--title \"$customRunName\"" : ''
     rfilename = customRunName ? "--filename " + customRunName.replaceAll('\\W','_').replaceAll('_+','_') + "_" + prefix + "_multiqc_report" : ''
     metadataOpts = params.metadata ? "--metadata ${metadata}" : ""
     splanOpts = params.samplePlan ? "--splan ${params.samplePlan}" : ""
     """
     awk -F"," -v sname=${prefix} '\$1==sname{print}' ${splan} > splan_${prefix}.csv
     stats2multiqc.sh splan_${prefix}.csv	
     mqc_header.py --name "${prefix} - nf-VIF" --version ${workflow.manifest.version} ${metadataOpts} ${splanOpts} > multiqc-config-header.yaml
     multiqc . -f $rtitle -n ${prefix}_nfvif_report.html -c $multiqcConfig -c 'mqc/hpv_config.mqc' -c multiqc-config-header.yaml -m fastqc -m custom_content
     """
   }
}else{
   process makeHpvConfig {

      when:
      !params.skipMultiqc

      input:
      file(geno) from hpvGenoMqcConfig.collect()
      file genes from chHpvGenesCoord.collect()

      output:
      file('*conf.mqc') into mqcHpvConf
      file('*bkp.mqc') into mqcGenepos

      script:
      prefix = geno.toString() - "_HPVgenotyping.filtered"
      """
      awk -F, '{print \$2}' ${geno} | sort -u > allgenotypes_unique.txt
      scrape_mqc_config.py  allgenotypes_unique.txt > hpv_conf.mqc
      gene_tracks.sh allgenotypes_unique.txt ${genes} 'genes' 
      """
   }

   process multiqcAllsamples {
     publishDir "${params.outdir}/MultiQC/", mode: 'copy'

     when:
     !params.skipMultiqc

     input:
     file splan from chSplan.first()
     file multiqcConfig from chMultiqcConfig.first()
     file hpvConfig from mqcHpvConf.collect().ifEmpty([])
     file('fastqc/*') from fastqcResults.map{items->items[1]}.collect().ifEmpty([])
     file('trimming/*') from trimgaloreResults.map{items->items[1]}.collect().ifEmpty([])
     file ('hpv/*') from hpvGenoStats.map{items->items[1]}.collect()
     file ('hpv/*') from hpvGenoMqc.map{items->items[1]}.collect()
     file ('hpv/*') from hpvCovStats.map{items->items[1]}.collect()
     file ('hpv/*') from hpvBwCov.map{items->items[1]}.collect()
     file ('hpv/*') from bkpPos.map{items->items[1]}.collect()
     file ('hpv/*') from mqcGenepos.map{items->items[1]}.collect()
     file ('hpv/*') from ttd.map{items->items[1]}.collect().ifEmpty([])
     file ('hpv/*') from hpvBowtie2Log.map{items->items[1]}.collect().ifEmpty([])
     file ('ctrl/*') from ctrlStats.map{items->items[1]}.collect()
  
     file ('software_versions/*') from software_versions_yaml.collect()
     file ('workflow_summary/*') from workflow_summary_yaml.collect()
 
     output:
     file splan
     file "*multiqc_report.html" into multiqcReport
     file "*_data"

     script:
     rtitle = customRunName ? "--title \"$customRunName\"" : ''
     rfilename = customRunName ? "--filename " + customRunName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
     metadataOpts = params.metadata ? "--metadata ${metadata}" : ""
     splanOpts = params.samplePlan ? "--splan ${params.samplePlan}" : ""
     """	
     stats2multiqc.sh ${splan}
     mqc_header.py --name "nf-VIF" --version ${workflow.manifest.version} ${metadataOpts} ${splanOpts} > multiqc-config-header.yaml         
     multiqc . -f $rtitle $rfilename -c $multiqcConfig -c $hpvConfig -c multiqc-config-header.yaml -m fastqc -m custom_content
     """
   }
}


