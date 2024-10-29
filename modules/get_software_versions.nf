process GET_SOFTWARE_VERSIONS {

    container = "phinguyen2000/hpv_version:9c95c92"

    input:
    path(scrape_software_versions_script)

    output:
    path('software_versions_mqc.yaml')
    // file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    bowtie2 --version > v_bowtie2.txt
    samtools --version > v_samtools.txt
    bedtools --version > v_bedtools.txt
    deeptools --version 2> v_deeptools.txt
    echo "BLAT v. 35" > v_blat.txt
    python --version 2> v_python.txt
    multiqc --version > v_multiqc.txt
    python $scrape_software_versions_script > software_versions_mqc.yaml
    """
}