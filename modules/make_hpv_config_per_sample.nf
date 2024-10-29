process MAKE_HPV_CONFIG_PER_SAMPLE {        
    tag "$prefix"
    container "phinguyen2000/pandas:813ad74"

    when:
    !params.skipMultiqc

    input:
    tuple path(sel_hpv_geno), path(hpv_genes_coord), path(scrape_mqc_config_script), path(gene_tracks_script)
    // file(geno) from hpvGenoMqcConfig
    // file genes from chHpvGenesCoord.collect()

    output:
    tuple val(prefix), path('*conf.mqc'),                 emit: mqc_hpv_conf
    tuple val(prefix), path('*bkp.mqc'), optional: true,  emit: mqc_genepos
    // set val(prefix),file('*conf.mqc') into mqcHpvConf
    // set val(prefix),file('*bkp.mqc') optional true into mqcGenepos

    script:
    prefix = sel_hpv_geno.toString() - "_HPVgenotyping.filtered"
    """
    awk -F, '{print \$2}' ${sel_hpv_geno} | sort -u > allgenotypes_unique.txt
    python $scrape_mqc_config_script  allgenotypes_unique.txt > ${prefix}_hpv_conf.mqc
    bash $gene_tracks_script allgenotypes_unique.txt ${hpv_genes_coord} ${prefix} 
    """
}