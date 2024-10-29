include {  HPV_LOCAL_MAPPING        }             from            "../modules/hpv_local_mapping.nf"
include {  HPV_LOCAL_MAPPING_STATS  }             from            "../modules/hpv_local_mapping_stats.nf"
include {  HPV_COVERAGE             }             from            "../modules/hpv_coverage.nf"


workflow LOCAL_MAPPING {
    take:
    trim_fastq                                  // ([val(prefix), listpath(trimmed_fastq)])
    filtered_sel_hpv_geno                       // ([path(sel_hpv_geno)])
    bwt2_index_hpv_split                        // ([path(bwt2_index_split)])

    main:
    HPV_LOCAL_MAPPING(
        filtered_sel_hpv_geno.splitCsv(header: ["sample", "hpv"])
                             .map{[ it["sample"], it["hpv"] ]}      // ([val(prefix), val(hpv)])
                             .combine(trim_fastq, by: 0)            // ([val(prefix), val(hpv), listpath(trimmed_fastq)])
                             .combine(bwt2_index_hpv_split)         // ([val(prefix), val(hpv), listpath(trimmed_fastq), path(bwt2_index_hpv_split)])
    )
    HPV_LOCAL_MAPPING_STATS(
        HPV_LOCAL_MAPPING.out                   // (val(prefix), val(hpv), path(local_bam))
    )
    HPV_COVERAGE(
        HPV_LOCAL_MAPPING.out                   // (val(prefix), val(hpv), path(local_bam))
    )

    emit:
    hpv_soft_bam      = HPV_LOCAL_MAPPING.out               // [(val(prefix), val(hpv), path(local_bam))]
    hpv_cov_stats     = HPV_LOCAL_MAPPING_STATS.out         // [(val(prefix), path(coverage_stats))]
    hpv_bw_cov        = HPV_COVERAGE.out.hpv_bw_cov         // [(val(prefix), path(covmatrix_mqc))]

}