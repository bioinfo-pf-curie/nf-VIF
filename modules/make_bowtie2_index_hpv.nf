process MAKE_BOWTIE2_INDEX_HPV {
    
    container "phinguyen2000/hpv_version:9c95c92"

    publishDir path: { params.saveReference ? "${params.outdir}/references" : params.outdir },
              saveAs: { params.saveReference ? it : null }, mode: 'copy'

    input:
    path(hpv_fasta_for_index)
    // file fasta from hpvFastaForIndex

    output:
    path("bowtie2IndexHpv")       , emit: main
    path("bowtie2IndexHpvSplit")  , emit: split
    // file "bowtie2IndexHpv" into bwt2IndexHpv
    // file "bowtie2IndexHpvSplit" into bwt2IndexHpvSplit

    script:
    hpv_bwt2_base = hpv_fasta_for_index.toString() - ~/(\.fa)?(\.fasta)?(\.fas)?$/
    """
    mkdir bowtie2IndexHpv
    bowtie2-build ${hpv_fasta_for_index} bowtie2IndexHpv/${hpv_bwt2_base}

    mkdir fastaSplit && cd fastaSplit
    cat ../$hpv_fasta_for_index | awk '{ if (substr(\$0, 1, 1)==">") {filename=(substr(\$0,2) ".fa")} print \$0 > filename }'
    cd .. && ls fastaSplit/* > listoffasta.txt

    mkdir bowtie2IndexHpvSplit
    while read ff; do
    base="\$(basename \"\${ff}\" | sed -e 's/.fa//')"
    bowtie2-build "\${ff}" bowtie2IndexHpvSplit/"\${base}"
    done < listoffasta.txt
    """
}