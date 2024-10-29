include {  GET_SOFTWARE_VERSIONS       }                 from                "../modules/get_software_versions.nf"
include {  WORKFLOW_SUMMARY_MQC        }                 from                "../modules/workflow_summary_mqc.nf"
include {  MAKE_HPV_CONFIG_PER_SAMPLE  }                 from                "../modules/make_hpv_config_per_sample.nf"
include {  MULTIQC                     }                 from                "../modules/multiqc.nf"
include {  MAKE_HPV_CONFIG             }                 from                "../modules/make_hpv_config.nf"
include {  MULTIQC_ALL_SAMPLES         }                 from                "../modules/multiqc_all_samples.nf"


workflow MULTIQC_PROCESSING {
    take:
    vif_ob
    filtered_sel_hpv_geno                    // ([path(sel_hpv_geno)])
    ch_hpv_genes_coord
    ch_splan
    ch_multiqc_config
    fastqc_results                           // ([val(sname), path(fastqc_rs)])
    trimming_report                          // ([val(sname), path(trimming_report)])
    ctrl_stats                               // ([val(prefix), path(ctrl_stats)])
    hpv_bowtie2_log                          // ([val(prefix), path(hpv_bowtie2_log)])
    hpv_geno_mqc                             // ([val(prefix), path(sel_hpv_geno)])
    hpv_geno_stats                           // ([val(prefix), path(hpv_geno_stats)])
    hpv_cov_stats                            // [(val(prefix), path(coverage_stats))]
    hpv_bw_cov                               // [(val(prefix), path(covmatrix_mqc))]
    bkp_pos                                  // [(val(prefix), path(3prime_bkp_mqc)]
    ttd                                      // ([val(prefix), path(bkptable_filtered)])

    main:
    scrape_software_versions_script = Channel.fromPath("$projectDir/src/scrape_software_versions.py")
    GET_SOFTWARE_VERSIONS(scrape_software_versions_script)
    
    mqc_yaml_content = params.skipMultiqc? Channel.empty() : Channel.of(vif_ob.ymlContent)
    WORKFLOW_SUMMARY_MQC(mqc_yaml_content)


    scrape_mqc_config_script = Channel.fromPath("$projectDir/src/scrape_mqc_config.py")
    gene_tracks_script       = Channel.fromPath("$projectDir/src/gene_tracks.sh")
    stats2_multiqc_script    = Channel.fromPath("$projectDir/src/stats2multiqc.sh")
    mqc_header_script        = Channel.fromPath("$projectDir/src/mqc_header.py")




    MAKE_HPV_CONFIG_PER_SAMPLE(
        filtered_sel_hpv_geno.combine(ch_hpv_genes_coord)
                             .combine(scrape_mqc_config_script)
                             .combine(gene_tracks_script)
    )

    MAKE_HPV_CONFIG_PER_SAMPLE.out.mqc_hpv_conf
                                .join(fastqc_results, remainder: true)
                                .join(trimming_report, remainder: true)
                                .join(ctrl_stats)
                                .join(hpv_bowtie2_log)
                                .join(hpv_geno_mqc)
                                .join(hpv_cov_stats.groupTuple(), remainder: true)
                                .join(hpv_bw_cov.groupTuple(), remainder: true)
                                .join(bkp_pos.groupTuple(), remainder: true)
                                .join(hpv_geno_stats, remainder: true)
                                .join(MAKE_HPV_CONFIG_PER_SAMPLE.out.mqc_genepos, remainder: true)
                                .join(ttd.groupTuple(), remainder: true)
                                .dump(tag: "join")
                                .set{ch_hpv_report}

    if (params.splitReport){

    MULTIQC(
        ch_splan.first().combine(ch_multiqc_config.first())
                        .combine(ch_hpv_report)
                        .combine(GET_SOFTWARE_VERSIONS.out)
                        .combine(WORKFLOW_SUMMARY_MQC.out)
                        .combine(stats2_multiqc_script)
                        .combine(mqc_header_script)
    )
    }else{
    
    MAKE_HPV_CONFIG(
        filtered_sel_hpv_geno.combine(ch_hpv_genes_coord)
                             .combine(scrape_mqc_config_script)
                             .combine(gene_tracks_script)
    )
    qc_all_sample_ch = ch_splan.combine(ch_multiqc_config)
    MULTIQC_ALL_SAMPLES(
        ch_splan.first().combine(ch_multiqc_config.first().map{[it]})
                       .combine(MAKE_HPV_CONFIG_PER_SAMPLE.out.mqc_hpv_conf.map{items->items[1]}.collect().ifEmpty([]).map{[it]})
                       .combine(fastqc_results.map{items->items[1]}.collect().ifEmpty([]).map{[it]})
                       .combine(trimming_report.map{items->items[1]}.collect().ifEmpty([]).map{[it]})
                       .combine(hpv_geno_stats.map{items->items[1]}.collect().map{[it]})
                       .combine(hpv_geno_mqc.map{items->items[1]}.collect().map{[it]})
                       .combine(hpv_cov_stats.map{items->items[1]}.collect().map{[it]})
                       .combine(hpv_bw_cov.map{items->items[1]}.collect().map{[it]})
                       .combine(bkp_pos.map{items->items[1]}.collect().map{[it]})
                       .combine(MAKE_HPV_CONFIG_PER_SAMPLE.out.mqc_genepos.map{items->items[1]}.collect().map{[it]})
                       .combine(ttd.map{items->items[1]}.collect().ifEmpty([]).map{[it]})
                       .combine(hpv_bowtie2_log.map{items->items[1]}.collect().ifEmpty([]).map{[it]})
                       .combine(ctrl_stats.map{items->items[1]}.collect().map{[it]})
                       .combine(GET_SOFTWARE_VERSIONS.out.collect().map{[it]})
                       .combine(WORKFLOW_SUMMARY_MQC.out.collect().map{[it]})
                       .combine(stats2_multiqc_script)
                       .combine(mqc_header_script)

    )
    }


}