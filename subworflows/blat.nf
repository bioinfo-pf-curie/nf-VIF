include    {  EXTRACT_BREAKPOINTS_SEQUENCE  }               from        "../modules/extract_breakpoints_sequence.nf"
include    {  BUILD_TEST_BLAT_DB            }               from        "../modules/build_test_blat_db.nf"
include    {  BLAT_SOFT_CLIPPED_SEQ         }               from        "../modules/blat_soft_clipped_seq.nf"
include    {  BLAT_SUMMARY                  }               from        "../modules/blat_summary.nf"


workflow BLAT {
    take:
    reference_fasta_for_index
    hpv_soft_bam            // [(val(prefix), val(hpv), path(local_bam))]

    main:

    /*
    * Breakpoint detection
    */
    extract_softclipped_code = Channel.fromPath("$projectDir/src/extractSoftclipped.py")
    EXTRACT_BREAKPOINTS_SEQUENCE(
        hpv_soft_bam.combine(extract_softclipped_code)          // [(val(prefix), val(hpv), path(local_bam), path(extract_softclipped_code))]
    )

    if (!params.skipBlat && params.blatdb){
        Channel.fromPath( params.blatdb )
                .ifEmpty { exit 1, "BLAT database not found: ${params.blatdb}" }
                .set { blatDatabase }
    } 
    
    /*
    * BLAT
    */
    if (!params.skipBlat && !params.blatdb && (workflow.profile =~ /test/)){
        BUILD_TEST_BLAT_DB(reference_fasta_for_index)
        blatDatabase = BUILD_TEST_BLAT_DB.out
    }

    if(params.skipBlat){
        blatDatabase = Channel.empty()
    }


    BLAT_SOFT_CLIPPED_SEQ(
        EXTRACT_BREAKPOINTS_SEQUENCE.out.clipped_seq                                    // [(val(pfix), val(prefix), path(clipped_seq))]
                                    .combine(blatDatabase)                              // [(val(pfix), val(prefix), path(clipped_seq), path(blatdb))]
    )

    blat_parser_script = Channel.fromPath("$projectDir/src/blatParser.py")
    BLAT_SUMMARY(
        BLAT_SOFT_CLIPPED_SEQ.out
                        .combine(
                            EXTRACT_BREAKPOINTS_SEQUENCE.out.bkp_info,                  // [(val(pfix), path(bkp_info))]
                            by: 0)                              
                        .combine(blat_parser_script)            // [(val(pfix), val(prefix), path(blat_res), path(bkp_info), path(blat_parser_script)]
    )

    emit:
    ttd                 =   BLAT_SUMMARY.out.ttd                            // ([val(prefix), path(bkptable_filtered)])
    bkp_pos             =   EXTRACT_BREAKPOINTS_SEQUENCE.out.bkp_pos        // [(val(prefix), path(3prime_bkp_mqc)]
}