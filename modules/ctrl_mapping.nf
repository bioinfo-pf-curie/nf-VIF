process CTRL_MAPPING {
  tag "$prefix"

  publishDir "${params.outdir}/ctrlMapping/", mode: 'copy',
      saveAs: {filename ->
          if (filename.endsWith(".log")) "logs/$filename"
          else if (params.saveAlignedIntermediates) filename
	  else null}

  input:
  tuple val(prefix), path(trimmed_fastq), path(bwt2_index_ctrl_folder)
  //   set val(prefix), file(reads) from readsCtrl
  //   file index from  bwt2IndexCtrl.collect()

  output:
  tuple val(prefix), file("${prefix}_ctrl.bam")                 , emit: ctrl_bam
  tuple val(prefix), file("*ctrl_bowtie2.log")                  , emit: ctrl_bowtie2_log
//   set val(prefix), file("${prefix}_ctrl.bam") into ctrlBam
//   set val(prefix), file("*ctrl_bowtie2.log") into ctrlBowtie2Log

  script:
  if ( params.singleEnd ){
  """
  bowtie2 --rg-id BMG --rg SM:${prefix} \\
          --very-sensitive \\
          -p ${task.cpus} \\
          -x ${bwt2_index_ctrl_folder}/ctrlRegions \\
          -U ${trimmed_fastq} > ${prefix}_ctrl.bam 2> ${prefix}_ctrl_bowtie2.log
  """
  }else{
  """
  bowtie2 --rg-id BMG --rg SM:${prefix} \\
          --very-sensitive \\
          -p ${task.cpus} \\
          -x ${bwt2_index_ctrl_folder}/ctrlRegions \\
          -1 ${trimmed_fastq[0]} -2 ${trimmed_fastq[1]} > ${prefix}_ctrl.bam 2> ${prefix}_ctrl_bowtie2.log
  """
  }
}