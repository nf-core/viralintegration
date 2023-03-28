/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowViralintegration.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta, params.viral_fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
ch_igvjs_VIF             = file("$projectDir/assets/igvjs_VIF.html", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CAT_FASTA } from '../modules/local/cat_fasta'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { HOST        } from '../subworkflows/local/host'
include { PLUS        } from '../subworkflows/local/plus'
include { INSERTION_SITES } from '../subworkflows/local/insertion_sites'
include { VALIDATE } from '../subworkflows/local/validate'
include { CHIMERIC_READS } from '../subworkflows/local/chimeric_reads'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { METAPHLAN3_METAPHLAN3       } from '../modules/nf-core/metaphlan3/metaphlan3/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow VIRALINTEGRATION {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // SUBWORKFLOW: Align input reads against host genome.
    //

    HOST (
        INPUT_CHECK.out.reads,
        params.fasta,
        params.gtf
    )
    ch_versions = ch_versions.mix(HOST.out.versions)

    //
    // MODULE: Run MetaPhlAn3
    //

    METAPHLAN3_METAPHLAN3 (
        HOST.out.polya_trimmed,
        params.metaphlan_db
    )
    ch_versions = ch_versions.mix(METAPHLAN3_METAPHLAN3.out.versions)

    //
    // MODULE: CAT viral and host reference fastas.
    //

    CAT_FASTA (
        params.fasta,
        params.viral_fasta
    )
    ch_versions = ch_versions.mix(CAT_FASTA.out.versions)

    //
    // SUBWORKFLOW: Align filtered reads against combined host and viral reference.
    //

    PLUS (
        HOST.out.polya_trimmed,
        CAT_FASTA.out.plus_fasta,
        params.gtf
    )
    ch_versions = ch_versions.mix(PLUS.out.versions)

    //
    // SUBWORKFLOW: Map potential virus insertion sites and create virus infection evidence viewer.
    //

    INSERTION_SITES (
        PLUS.out.bam_bai_junction,
        params.fasta,
        params.viral_fasta,
        PLUS.out.sam_bam,
        PLUS.out.bai,
        ch_igvjs_VIF
    )
    ch_versions = ch_versions.mix(INSERTION_SITES.out.versions)

    //
    // SUBWORKFLOW: Validate potential chimeric reads.
    //

    VALIDATE (
        HOST.out.fastq,
        INSERTION_SITES.out.fasta_extract,
        PLUS.out.index
    )
    ch_versions = ch_versions.mix(VALIDATE.out.versions)

    //
    // SUBWORKFLOW: Refine chimeric read counts and create final summary report.
    //

    CHIMERIC_READS (
        VALIDATE.out.bam_bai,
        INSERTION_SITES.out.gtf_extract,
        INSERTION_SITES.out.fasta_extract,
        INSERTION_SITES.out.filtered_abridged,
        INSERTION_SITES.out.genome_abundance_plot,
        INSERTION_SITES.out.read_counts_image,
        INSERTION_SITES.out.read_counts_log_image,
        params.gtf,
        ch_igvjs_VIF
    )
    ch_versions = ch_versions.mix(CHIMERIC_READS.out.versions)

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowViralintegration.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowViralintegration.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
        .mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        .mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
        .mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
        .mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
        .mix(HOST.out.log_final.collect{it[1]}.ifEmpty([]))
        .mix(HOST.out.mqc_log.collect{it[1]}.ifEmpty([]))
        .mix(PLUS.out.log_final.collect{it[1]}.ifEmpty([]))
        .mix(VALIDATE.out.log_final.collect{it[1]}.ifEmpty([]))


    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)

    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
