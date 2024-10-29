include {  CTRL_MAPPING      }              from            "../modules/ctrl_mapping.nf"
include {  CTRL_STATS        }              from            "../modules/ctrl_stats.nf"
include {  HPV_MAPPING       }              from            "../modules/hpv_mapping.nf"
include {  SELECT_GENOTYPES  }              from            "../modules/select_genotypes.nf"

workflow MAPPING {
    take:
    trim_fastq                          // ([val(name), listpath(trimmed_fastq)])
    bwt2_index_ctrl                     // ([path(bwt2_index_ctrl_folder)])
    bwt2_index_hpv                      // ([path(bwt2_index)])
    hpv_bwt2_base

    main:
    /*
    * Mapping on control regions
    */

    CTRL_MAPPING(
        trim_fastq.combine(bwt2_index_ctrl)        // ([val(prefix), listpath(trimmed_fastq), path(bwt2_index_ctrl_folder)])
    )
    CTRL_STATS(
        CTRL_MAPPING.out.ctrl_bam                  // ([val(prefix), path(ctrl_bam)])
    )

    /*
    * HPV mapping and genotyping
    */ 

    HPV_MAPPING(
        trim_fastq.combine(bwt2_index_hpv)          // ([val(prefix), listpath(trimmed_fastq), path(bwt2_index_hpv)])
                  .combine(hpv_bwt2_base)           // ([val(prefix), listpath(trimmed_fastq), path(bwt2_index_hpv), val(hpv_bwt2_base)])
    )
    SELECT_GENOTYPES(
        HPV_MAPPING.out.hpv_bam                     // ([val(prefix), path(hpvs_bam)])
    )

    emit:
    sel_hpv_geno     = SELECT_GENOTYPES.out.sel_hpv_geno
    ctrl_stats       = CTRL_STATS.out.ctrl_stats                // ([val(prefix), path(ctrl_stats)])
    hpv_bowtie2_log  = HPV_MAPPING.out.hpv_bowtie2_log          // ([val(prefix), path(hpv_bowtie2_log)])
    hpv_geno_mqc     = SELECT_GENOTYPES.out.hpv_geno_mqc        // ([val(prefix), path(sel_hpv_geno)])
    hpv_geno_stats   = SELECT_GENOTYPES.out.hpv_geno_stats      // ([val(prefix), path(hpv_geno_stats)])
}