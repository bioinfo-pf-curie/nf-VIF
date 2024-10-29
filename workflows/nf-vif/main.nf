include {  PREPROCESSING                 }               from        "../../subworflows/preprocessing.nf"
include {  QC                            }               from        "../../subworflows/qc.nf"
include {  MAPPING                       }               from        "../../subworflows/mapping.nf"
include {  LOCAL_MAPPING                 }               from        "../../subworflows/local_mapping.nf"
include {  BLAT                          }               from        "../../subworflows/blat.nf"
include {  MULTIQC_PROCESSING            }               from        "../../subworflows/multiqc_processing.nf"

workflow NF_VIF{
    take:
    bwt2_ref_index
    reference_fasta_for_index
    hpv_fasta_for_index
    ch_fasta_ctrl
    reads_trimgalore                        // ([val(prefix), listpath(fastq_file)])
    hpv_bwt2_base
    vif_ob
    ch_hpv_genes_coord
    ch_splan
    ch_multiqc_config


    main:
    /****************************************************
    * PRE-PROCESSING
    */

    PREPROCESSING(
        bwt2_ref_index,
        reference_fasta_for_index,
        hpv_fasta_for_index,
        ch_fasta_ctrl
    )

    /**************************** Main worflow ******************************


    /*
    * Quality control
    */
    QC(reads_trimgalore)

    /*
    * Mapping
    */

    MAPPING(
        QC.out.trim_fastq,                          // ([val(prefix), listpath(trimmed_fastq)])
        PREPROCESSING.out.bwt2_index_ctrl,          // ([path(bwt2_index_ctrl_folder)])
        PREPROCESSING.out.bwt2_index_hpv,           // ([path(bwt2_index)])
        hpv_bwt2_base
    )

    filtered_sel_hpv_geno = MAPPING.out.sel_hpv_geno
                                   .filter { geno -> vif_ob.checkGenotypes(geno) }                  // ([path(sel_hpv_geno)])
            // .into { hpvGenoFilter; hpvGenoMqcConfig }

    /*
    * Local Mapping for selected genotypes
    */

    LOCAL_MAPPING(
        QC.out.trim_fastq,                          // ([val(prefix), listpath(trimmed_fastq)])
        filtered_sel_hpv_geno,                      // ([path(sel_hpv_geno)])
        PREPROCESSING.out.bwt2_index_hpv_split      // ([path(bwt2_index_split)])
    )

    /*
    * Blat
    */  

    BLAT(
        reference_fasta_for_index,
        LOCAL_MAPPING.out.hpv_soft_bam,                          // [(val(prefix), val(hpv), path(local_bam))]
    )
    // ttd = Channel.from(false)


    /*
    /* MultiQC PROCESSING
    */

    MULTIQC_PROCESSING(
        vif_ob,
        filtered_sel_hpv_geno,
        ch_hpv_genes_coord,
        ch_splan,
        ch_multiqc_config,
        QC.out.fastqc_results,                       // ([val(sname), path(fastqc_rs)])
        QC.out.trimming_report,                      // ([val(sname), path(trimming_report)])
        MAPPING.out.ctrl_stats,                      // ([val(prefix), path(ctrl_stats)])
        MAPPING.out.hpv_bowtie2_log,                 // ([val(prefix), path(hpv_bowtie2_log)])
        MAPPING.out.hpv_geno_mqc,                    // ([val(prefix), path(sel_hpv_geno)])
        MAPPING.out.hpv_geno_stats,                  // ([val(prefix), path(hpv_geno_stats)])
        LOCAL_MAPPING.out.hpv_cov_stats,             // [(val(prefix), path(coverage_stats))]
        LOCAL_MAPPING.out.hpv_bw_cov,                // [(val(prefix), path(covmatrix_mqc))]
        BLAT.out.bkp_pos,                            // [(val(prefix), path(3prime_bkp_mqc)]
        BLAT.out.ttd                                 // ([val(prefix), path(bkptable_filtered)])
    )

}