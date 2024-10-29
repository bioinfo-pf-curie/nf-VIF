include {  MAKE_BOWTIE2_INDEX      }        from        "../modules/make_bowtie2_index.nf"
include {  MAKE_BOWTIE2_INDEX_HPV  }        from        "../modules/make_bowtie2_index_hpv.nf"
include {  MAKE_BOWTIE2_INDEX_CTRL }        from        "../modules/make_bowtie2_index_ctrl.nf"


workflow PREPROCESSING{
    take:
    bwt2_ref_index
    reference_fasta_for_index
    hpv_fasta_for_index
    ch_fasta_ctrl

    main:
    MAKE_BOWTIE2_INDEX(reference_fasta_for_index)
    MAKE_BOWTIE2_INDEX_HPV(hpv_fasta_for_index)
    MAKE_BOWTIE2_INDEX_CTRL(ch_fasta_ctrl)

    emit:
    bwt2_index_ctrl      = MAKE_BOWTIE2_INDEX_CTRL.out                  // ([path(bwt2_index_ctrl_folder)])
    bwt2_index_hpv       = MAKE_BOWTIE2_INDEX_HPV.out.main              // ([path(bwt2_index)])
    bwt2_index_hpv_split = MAKE_BOWTIE2_INDEX_HPV.out.split             // ([path(bwt2_index_split)])
}

