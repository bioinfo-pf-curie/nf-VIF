process BLAT_SUMMARY {
    tag "$prefix:$pfix"
    container "phinguyen2000/pandas:b4381b3"

    publishDir "${params.outdir}/hpvMapping/blat", mode: 'copy'

    input:
    tuple val(pfix), val(prefix), path(blat_res), path(bkp_info), path(blat_parser_script)
    // set val(pfix), val(sname), file(psl), file(csv) from blatRes.join(bkpInfo).dump(tag:"blat")

    output:
    tuple val(prefix), path("*_table_filtered.csv"),                 emit:  ttd_hq
    tuple val(prefix), path("*_table.csv"),                          emit:  table_hq
    tuple val(prefix), path("*_bkptable_filtered.csv"),              emit:  ttd
    tuple val(prefix), path("*_bkptable.csv"),                       emit:  table


    // set val(prefix), file("*_table_filtered.csv") into ttdHQ
    // set val(prefix), file("*_table.csv") into tableHQ
    // set val(prefix), file("*_bkptable_filtered.csv") into ttd
    // set val(prefix), file("*_bkptable.csv") into table

    script:
    """
    python $blat_parser_script -f ${blat_res} -b ${bkp_info} --sname ${prefix}
    """
}