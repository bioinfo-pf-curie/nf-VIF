process MAKE_HPV_CONFIG {
    tag "$prefix"
    container = "phinguyen2000/pandas:813ad74"

    when:
    !params.skipMultiqc

    input:
    tuple path(sel_hpv_geno), path(hpv_genes_coord), path(scrape_mqc_config_script), path(gene_tracks_script)
    // file(geno) from hpvGenoMqcConfig.collect()
    // file genes from chHpvGenesCoord.collect()

    output:
    path('*conf.mqc'),                  emit: mqc_hpv_conf
    path('*bkp.mqc'),                   emit: mqc_genepos
    // file('*conf.mqc') into mqcHpvConf
    // file('*bkp.mqc') into mqcGenepos

    script:
    prefix = sel_hpv_geno.toString() - "_HPVgenotyping.filtered"
    """
    awk -F, '{print \$2}' ${sel_hpv_geno} | sort -u > allgenotypes_unique.txt
    python $scrape_mqc_config_script  allgenotypes_unique.txt > hpv_conf.mqc
    bash $gene_tracks_script allgenotypes_unique.txt ${hpv_genes_coord} 'genes' 
    """
}