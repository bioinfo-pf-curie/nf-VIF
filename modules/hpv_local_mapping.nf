process HPV_LOCAL_MAPPING {
  tag "$prefix:$hpv"

  publishDir "${params.outdir}/hpvMapping/pergenotype", mode: 'copy',
      saveAs: {filename ->
          if (filename.endsWith(".log")) "logs/$filename"
	  else if (params.saveAlignedIntermediates) filename
          else null}
  
  input:
  tuple val(prefix), val(hpv), path(trimmed_fastq), path(bwt2_index_hpv_split)

  // file index from bwt2IndexHpvSplit.first()
  // set val(prefix), val(hpv), file(reads) from hpvGenoFilter
  //   .splitCsv(header: ["sample", "hpv"])
  //   .map{
  //     [ it["sample"], it["hpv"] ]
  //   }
  //   .combine(readsSplitmap, by: 0)
  //   .dump(tag: "hpvloc")

  output:
  tuple val(prefix), val(hpv), path("${prefix}-${hpv}.bam")
  // set val(prefix), file("*.bam") into hpvLocalBam, hpvCovBam, hpvSoftBam

  script: 
  if ( params.singleEnd ){
  """
  bowtie2 --rg-id BMG --rg SM:${prefix} \\
          --local --very-sensitive-local --no-unal \\
          -p ${task.cpus} \\
          -x ${bwt2_index_hpv_split}/${hpv} \\
          -U ${trimmed_fastq} > ${prefix}-${hpv}.bam 2> ${prefix}-${hpv}_bowtie2.log
  """
  }else{
  """
  bowtie2 --rg-id BMG --rg SM:${prefix} \\
          --local --very-sensitive-local --no-unal \\
          -p ${task.cpus} \\
          -x ${bwt2_index_hpv_split}/${hpv} \\
          -1 ${trimmed_fastq[0]} -2 ${trimmed_fastq[1]} > ${prefix}-${hpv}.bam 2> ${prefix}-${hpv}_bowtie2.log
  """
  }
}