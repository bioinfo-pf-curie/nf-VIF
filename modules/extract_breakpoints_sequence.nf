process EXTRACT_BREAKPOINTS_SEQUENCE {
   tag "$prefix:$hpv"
   container "phinguyen2000/pandas:b4381b3"

   publishDir "${params.outdir}/hpvMapping/softclipped", mode: 'copy',
               saveAs: {filename -> 
                   if (filename.indexOf(".mqc") > 0) "mqc/$filename"
		   else filename}

   input:
   tuple val(prefix), val(hpv), path(local_bam), path(extract_softclipped_code)
   // set val(prefix), file(bam) from hpvSoftBam

   output:
   tuple val(prefix), path("${pfix}_3prime_bkp.mqc"),            emit: bkp_pos
   tuple val(pfix), path("*.csv"),                               emit: bkp_info
   tuple val(pfix), val(prefix), path("*.fa"),                   emit: clipped_seq

   // set val(prefix), file("*.mqc") into bkpPos mode 'flatten'
   // set val(pfix), file("*.csv") into bkpInfo
   // set val(pfix), val(prefix), file("*.fa") into clippedSeq

   script:
   pfix= local_bam.toString() - ~/.bam$/
   """
   python $extract_softclipped_code -v --mqc --stranded --minLen ${params.minLen} ${local_bam}
   sort -k1,1n ${pfix}_3prime_bkp.mqc | awk 'BEGIN{nr=1} nr<\$1{for (i=nr;i<\$1;i++){print i"\t"0} nr=\$1}{print; nr+=1}' > file1.tmp
   mv file1.tmp ${pfix}_3prime_bkp.mqc
   sort -k1,1n ${pfix}_5prime_bkp.mqc | awk 'BEGIN{nr=1} nr<\$1{for (i=nr;i<\$1;i++){print i"\t"0} nr=\$1}{print; nr+=1}' > file2.tmp
   mv file2.tmp ${pfix}_5prime_bkp.mqc
   """
}