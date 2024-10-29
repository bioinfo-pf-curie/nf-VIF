process HPV_COVERAGE {
  tag "$prefix:$hpv"
  container "phinguyen2000/deeptools-samtools:c413f4f"
  
  input:
  tuple val(prefix), val(hpv), path(local_bam)
  // set val(prefix), file(bam) from hpvCovBam

  output:
  tuple val(prefix), path("${pfix}_covmatrix.mqc"),             emit: hpv_bw_cov
  tuple val(prefix), path("${pfix}_sorted.{bam,bam.bai}"),      emit: hpv_sorted_bams
  // set val(prefix), file("*covmatrix.mqc") into hpvBwCov
  // set val(prefix), file('*sorted.{bam,bam.bai}') into hpvSortedBams

  script:
  pfix= local_bam.toString() - ~/(_sorted)?(.bam)?$/
  normOpts = params.splitReport ? "" : "--normalizeUsing CPM"
  """
  samtools sort -@ ${task.cpus} -o ${pfix}_sorted.bam ${local_bam}
  samtools index ${pfix}_sorted.bam
  bamCoverage -b ${pfix}_sorted.bam --binSize 50 ${normOpts} --outFileFormat bedgraph -o ${pfix}.bedgraph
  awk -F"\t" '{OFS="\t"; print \$2+25,\$4}' ${pfix}.bedgraph > ${pfix}_covmatrix.mqc
  """
}