process MULTIQC {
    tag "$prefix"
    container "phinguyen2000/multiqc:c42a7c6"

    publishDir "${params.outdir}/MultiQC/", mode: 'copy'

    when:
    !params.skipMultiqc

    input:
    tuple path(splan), path(multiqc_config), 
            val(prefix), path('mqc/hpv_config.mqc'), path('fastqc/*'), path('trimming/*'), path('ctrl/*'),
            path('hpv/*'), path('hpv/*'), path('hpv/*'), path('hpv/*'), path('hpv/*'), path('hpv/*'), path('hpv/*'), path('hpv/*'),
            path('software_versions/*'),
            path('workflow_summary/*'),
            path(stats2_multiqc_script),
            path(mqc_header_script)

    // file splan from chSplan.first()
    // file multiqcConfig from chMultiqcConfig.first()
    // set val(prefix), file('mqc/hpv_config.mqc'), file('fastqc/*'), file('trimming/*'), file('ctrl/*'), 
    // file('hpv/*'), file('hpv/*'), file('hpv/*'), file('hpv/*'), file('hpv/*'), file('hpv/*'), file('hpv/*'), file('hpv/*') from chHpvReport.dump(tag: "mqc") 
    // file ('software_versions/*') from software_versions_yaml.collect()
    // file ('workflow_summary/*') from workflow_summary_yaml.collect()

    output:
    path "*.html",          emit: multiqc_report
    path "*_data",          emit: data

    // file splan
    // file "*.html" into multiqc_report
    // file "*_data"

    script:
    rtitle = workflow.runName ? "--title \"${workflow.runName}\"" : ''
    rfilename = workflow.runName ? "--filename " + workflow.runName.replaceAll('\\W','_').replaceAll('_+','_') + "_" + prefix + "_multiqc_report" : ''
    metadataOpts = params.metadata ? "--metadata ${metadata}" : ""
    splanOpts = params.samplePlan ? "--splan ${params.samplePlan}" : ""
    """
    awk -F"," -v sname=${prefix} '\$1==sname{print}' ${splan} > splan_${prefix}.csv
    bash $stats2_multiqc_script splan_${prefix}.csv	
    python $mqc_header_script --name "${prefix} - nf-VIF" --version ${workflow.manifest.version} ${metadataOpts} ${splanOpts} > multiqc-config-header.yaml
    multiqc . -f $rtitle -n ${prefix}_nfvif_report.html -c $multiqc_config -c 'mqc/hpv_config.mqc' -c multiqc-config-header.yaml -m fastqc -m custom_content
    """
}