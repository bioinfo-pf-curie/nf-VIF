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


/*
 * SET UP CONFIGURATION VARIABLES
 */

nextflow.enable.dsl = 2

def vif_ob = new VIF(workflow, log, params, "$USER")

// Show help emssage
if (params.help){
    vif_ob.helpMessage()
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
custom_run_name = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_run_name = workflow.runName
}

// Stage config files
ch_multiqc_config = Channel.fromPath(params.multiqc_config)
ch_output_docs = Channel.fromPath("$baseDir/docs/output.md")
ch_fasta_ctrl = Channel.fromPath(params.fastaCtrl)
ch_hpv_genes_coord = Channel.fromPath(params.genesHpv)

/*
 * CHANNELS
 */

/*
 * Create a channel for input read files
 */

if(params.samplePlan){
   if(params.singleEnd){
      Channel
         .fromPath("${params.samplePlan}")
         .splitCsv(header: false)
         .map{ row -> [ row[0], [row[2]]] }
         .set {reads_trimgalore}
   }else{
      Channel
         .fromPath("${params.samplePlan}")
         .splitCsv(header: false)
         .map{ row -> [ row[0], [row[2], row[3]]] }
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
       .collectFile{
         item -> ["samplePlan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + '\n']
        }
       .set{ ch_splan }
  }else{
     Channel
       .from(params.readPaths)
       .collectFile{
         item -> ["samplePlan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + ',' + item[1][1] + '\n']
        }
       .set{ ch_splan }
  }
}else{
  if (params.singleEnd){
    Channel
       .fromFilePairs( params.reads, size: 1 )
       .collectFile() {
          item -> ["samplePlan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + '\n']
       }     
       .set { ch_splan }
  }else{
    Channel
       .fromFilePairs( params.reads, size: 2 )
       .collectFile() {
          item -> ["samplePlan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + ',' + item[1][1] + '\n']
       }     
       .set { ch_splan }
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
      .set { bwt2_ref_index }
   
   reference_fasta_for_index = Channel.empty()
}
else if ( params.fasta ) {
   lastPath = params.fasta.lastIndexOf(File.separator)
   refBwt2Base = params.fasta.substring(lastPath+1)

   bwt2_ref_index = Channel.empty()
   Channel.fromPath( params.fasta )
        .ifEmpty { exit 1, "Genome index: Fasta file not found: ${params.fasta}" }
        .set { reference_fasta_for_index }
}
else {
   exit 1, "No reference genome specified!"
}

//HPV genome

if ( params.bwt2IndexHpv && params.bwt2IndexHpvSplit ){

   lastPath = params.bwt2IndexHpv.lastIndexOf(File.separator)
   hpvBwt2Dir =  params.bwt2IndexHpv.substring(0,lastPath+1)
   hpv_bwt2_base = Channel.of(params.bwt2IndexHpv.substring(lastPath+1))

   Channel.fromPath( hpvBwt2Dir , checkIfExists: true)
      .ifEmpty { exit 1, "HPV index: Provided index not found: ${params.bwt2IndexHpv}" }
      .set { bwt2IndexHpv }

   lastPath = params.bwt2IndexHpvSplit.lastIndexOf(File.separator)
   hpvSplitBwt2Dir =  params.bwt2IndexHpvSplit.substring(0,lastPath+1)
 
   Channel.fromPath( hpvSplitBwt2Dir , checkIfExists: true)
      .ifEmpty { exit 1, "HPV index per strain: Provided index not found: ${params.bwt2IndexHpvSplit}" }
      .set { bwt2IndexHpvSplit }

   hpv_fasta_for_index = Channel.empty()
}
else if ( params.fastaHpv ){
   lastPath = params.fastaHpv.lastIndexOf(File.separator)
   hpv_bwt2_base = Channel.of(params.fastaHpv.substring(lastPath+1) - ~/(\.fa)?(\.fasta)?(\.fas)?$/)

   Channel.fromPath( params.fastaHpv )
        .ifEmpty { exit 1, "HPV index: Fasta file not found: ${params.fastaHpv}" }
        .set { hpv_fasta_for_index }
}
else{
   exit 1, "No HPV genome specified!"
}

// Header log info
vif_ob.headerInfo()


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MAIN WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include {  NF_VIF  }                from           "./workflows/nf-vif/main.nf"

workflow{
   NF_VIF(
      bwt2_ref_index,
      reference_fasta_for_index,
      hpv_fasta_for_index,
      ch_fasta_ctrl,
      reads_trimgalore,                        // ([val(prefix), listpath(fastq_file)])
      hpv_bwt2_base,
      vif_ob,
      ch_hpv_genes_coord,
      ch_splan,
      ch_multiqc_config
   )
}


