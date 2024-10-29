process HPV_MAPPING {
  tag "$prefix"

  publishDir "${params.outdir}/hpvMapping/allref", mode: 'copy',
        saveAs: {filename ->
            if (filename.endsWith(".log")) "logs/$filename"
            else if (params.saveAlignedIntermediates) filename
	    else null}

  input:
  tuple val(prefix), path(trimmed_fastq), path(bwt2_index_hpv), val(hpv_bwt2_base)
//   set val(prefix), file(reads) from readsHpvmap
//   file index from bwt2IndexHpv.collect()

  output:
  tuple val(prefix), path("${prefix}_hpvs.bam"),                emit:  hpv_bam
  tuple val(prefix), path("*bowtie2.log"),                      emit:  hpv_bowtie2_log
//   set val(prefix), file("${prefix}_hpvs.bam") into hpvBam
//   set val(prefix), file("*bowtie2.log") into hpvBowtie2Log

  script:
  if ( params.singleEnd ){
  """
  bowtie2 --rg-id BMG --rg SM:${prefix} \\
          --very-sensitive --no-unal \\
          -p ${task.cpus} \\
          -x ${bwt2_index_hpv}/${hpv_bwt2_base} \\
          -U ${trimmed_fastq} > ${prefix}_hpvs.bam 2> ${prefix}_hpvs_bowtie2.log
  """
  }else{
  """
  bowtie2 --rg-id BMG --rg SM:${prefix} \\
          --very-sensitive --no-unal \\
          -p ${task.cpus} \\
          -x ${bwt2_index_hpv}/${hpv_bwt2_base} \\
          -1 ${trimmed_fastq[0]} -2 ${trimmed_fastq[1]} > ${prefix}_hpvs.bam 2> ${prefix}_hpvs_bowtie2.log 
  """
  }
}