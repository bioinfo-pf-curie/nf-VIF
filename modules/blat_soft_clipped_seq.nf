process BLAT_SOFT_CLIPPED_SEQ {
    tag "$prefix:$pfix"

    publishDir "${params.outdir}/hpvMapping/blat", mode: 'copy'

    container "phinguyen2000/blat:255336f"

    input:
    tuple val(pfix), val(prefix), path(clipped_seq), path(blatdb)
    // file(blatdb) from blatDatabase.collect()
    // set val(pfix), val(sname), file(fasta) from clippedSeq

    output:
    tuple val(pfix), val(prefix), path("*.tsv")
    // set val(pfix), val(sname), file("*.tsv") into blatRes

    script:
    """
    blat ${blatdb} ${clipped_seq} ${pfix}.tsv -noHead -minScore=25 -minIdentity=90
    """
}