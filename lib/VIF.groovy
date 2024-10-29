class VIF {

    def workflow
    def log
    def params
    def user
    Map summary = [:]
    public String ymlContent

    VIF(workflow, log, params, user){
        this.workflow = workflow
        this.log = log
        this.params = params
        this.user = user
        initializeSummary()
        this.ymlContent = initializeYamlContent()
    }


    public void helpMessage() {
        this.log.info"""
        
        HPV v${this.workflow.manifest.version}
        =======================================================

        Usage:

        nextflow run main.nf --reads '*_R{1,2}.fastq.gz' --genome 'hg19' -profile conda
        nextflow run main.nf --samplePlan sample_plan --genome hg19 -profile conda

        Mandatory arguments:
        --reads                       Path to input data (must be surrounded with quotes)
        --samplePlan                  Path to sample plan file if '--reads' is not specified
        --genome                      Name of iGenomes reference
        -profile                      Configuration profile to use. test / toolsPath / singularity / cluster (see below)

        Options:
        --singleEnd                   Specifies that the input is single end reads

        Genome References:              If not specified in the configuration file or you wish to overwrite any of the references.
        --genome                      Name of iGenomes reference
        --bwt2Index                   Path to Bowtie2 index
        --fasta                       Path to Fasta reference (.fasta)
        --blatdb                      Path to BLAT database (.2bit)

        HPV References:
        --fastaHpv                    Path to Fasta HPV reference (.fasta)                 
        --bwt2IndexHpv                Path to Bowtie2 index for all HPV strains
        --bwt2IndexHpvSplit           Path to Bowtie2 index per HPV strain
        --saveReference               Save all references generated during the analysis. Default: False

        Advanced options:
        --minMapq                     Minimum reads mapping quality. Default: 0
        --minLen                      Minimum trimmed length sequence to consider. Default: 15
        --minFreqGeno                 Fraction of reads to consider a genotpye. Default: 0.2
        --nbGeno                      Number of HPV genotype to consider
        --splitReport                 Generate one report per sample
    
        Other options:
        --outdir                      The output directory where the results will be saved
        --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
        -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

        Skip options:
        --skipTrimming               Skip trimming step
        --skipFastqc                 Skip quality controls on sequencing reads
        --skipBlat                   Skip Human mapping with Blat
        --skipMultiqc                Skip report

        =======================================================                                                                                                                                                 
        Available Profiles
        -profile test                Set up the test dataset
        -profile toolsPath           Use the paths defined in configuration for each tool
        -profile singularity         Use the Singularity images for each process
        -profile cluster             Run the workflow on the cluster, instead of locally
        -profile gcloud              Run the workflow with google Batch

        """.stripIndent()
    }

    
    private void initializeSummary() {
        summary['Pipeline Name']  = 'HPV'
        summary['Pipeline Version'] = this.workflow.manifest.version
        summary['Run Name']     = this.workflow.runName
        if (this.params.samplePlan) {
        summary['SamplePlan']   = this.params.samplePlan
        }else{
        summary['Reads']        = this.params.reads
        }
        summary['Fasta Ref']      = this.params.fasta
        summary['BLAT database']  = this.params.blatdb
        summary['Fasta HPV']      = this.params.fastaHpv
        summary['Min MAPQ']       = this.params.minMapq
        summary['Min length']     = this.params.minLen
        summary['Min Freq Geno']  = this.params.minFreqGeno
        summary['Split report']   = this.params.splitReport
        summary['Max Memory']     = this.params.max_memory
        summary['Max CPUs']       = this.params.max_cpus
        summary['Max Time']       = this.params.max_time
        summary['Output dir']     = this.params.outdir
        summary['Working dir']    = this.workflow.workDir
        summary['Container Engine'] = this.workflow.containerEngine
        summary['Current user']   = this.user
        summary['Config Profile'] = this.workflow.profile
        if(this.params.email) {
            summary['E-mail Address'] = this.params.email
        }
    }

    public void headerInfo(){
        this.log.info """=======================================================

HPV v${this.workflow.manifest.version}
======================================================="""

        
        this.log.info this.summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
        this.log.info "======================================================="
    }

    public boolean checkGenotypes(geno){
        // Filter - removes all samples for which the genotype has not been detected
        def skippedNogeno = []
        def nbGeno = 0;
        geno.eachLine { nbGeno++; }
        def samplename = geno.getBaseName() - '_HPVgenotyping.filered'
        if(nbGeno < 1 ){
            this.log.info "#################### NO HPV GENOTYPE DETECTED! IGNORING FOR FURTHER DOWNSTREAM ANALYSIS! ($samplename)"
            skippedNogeno << samplename
            return false
        } else {
            return true
        }
    }

    public String initializeYamlContent(){
        String context = """
        id: 'summary'
        description: \\" - this information is collected when the pipeline is started.\\"
        section_name: 'Workflow Summary'
        section_href: 'https://gitlab.curie.fr/illumina-hpv'
        plot_type: 'html'
        data: |
            <dl class=\\"dl-horizontal\\">
            ${this.summary.collect { k, v -> "                <dt>$k</dt><dd><samp>${v ?: '<span style=\\"color:#999999;\\">N/A</span>'}</samp></dd>" }.join("\n")}
            </dl>
        """.stripIndent()

        return context
    }
}