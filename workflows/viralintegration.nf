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

include { POLYA_STRIPPER } from '../modules/local/polyA_stripper'
include { CAT_FASTA } from '../modules/local/cat_fasta'
include { INSERTION_SITE_CANDIDATES } from '../modules/local/insertion_site_candidates'
include { ABRIDGED_TSV } from '../modules/local/abridged_tsv'
include { VIRUS_REPORT } from '../modules/local/virus_report'
include { EXTRACT_CHIMERIC_GENOMIC_TARGETS } from '../modules/local/extract_chimeric_genomic_targets'
include { STAR_ALIGN_VALIDATE } from '../modules/local/star_align_validate'
include { CHIMERIC_CONTIG_EVIDENCE_ANALYZER } from '../modules/local/chimeric_contig_evidence_analyzer'
include { SUMMARY_REPORT } from '../modules/local/summary_report'
include { REMOVE_DUPLICATES } from '../modules/local/remove_duplicates'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { HOST        } from '../subworkflows/local/host'
include { PLUS        } from '../subworkflows/local/plus'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_VALIDATE
          SAMTOOLS_SORT as SAMTOOLS_SORT_DUPLICATES } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_VALIDATE
          SAMTOOLS_INDEX as SAMTOOLS_INDEX_DUPLICATES } from '../modules/nf-core/samtools/index/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
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

    INSERTION_SITE_CANDIDATES (
        PLUS.out.bam_bai_junction,
        params.fasta,
        params.viral_fasta
    )
    ch_versions = ch_versions.mix(INSERTION_SITE_CANDIDATES.out.versions.first())

    ABRIDGED_TSV (
        INSERTION_SITE_CANDIDATES.out.full
    )
    // TODO ch_versions = ch_versions.mix(ABRIDGED_TSV.out.versions.first())

    PLUS.out.sam_bam
        .join(PLUS.out.bai, by: [0], remainder: true)
        .join(ABRIDGED_TSV.out.filtered_abridged)
        .set { ch_bam_bai_filtered }

    VIRUS_REPORT (
        ch_bam_bai_filtered,
        params.viral_fasta,
        ch_igvjs_VIF
    )
    ch_versions = ch_versions.mix(VIRUS_REPORT.out.versions.first())

    // TODO Handle insertion_site_candidates
    // File insertion_site_candidates_use = select_first([insertion_site_candidates, InsertionSiteCandidates.filtered_abridged])

    EXTRACT_CHIMERIC_GENOMIC_TARGETS (
        ABRIDGED_TSV.out.filtered_abridged,
        params.fasta,
        params.viral_fasta
    )
    ch_versions = ch_versions.mix(EXTRACT_CHIMERIC_GENOMIC_TARGETS.out.versions.first())

    HOST.out.fastq
        .join(EXTRACT_CHIMERIC_GENOMIC_TARGETS.out.fasta_extract)
        .set { ch_unaligned_fastq_fasta }
    STAR_ALIGN_VALIDATE (
        ch_unaligned_fastq_fasta,
        PLUS.out.index,
        "illumina",
        false
    )
    ch_versions = ch_versions.mix(STAR_ALIGN_VALIDATE.out.versions.first())

    SAMTOOLS_SORT_VALIDATE (
        STAR_ALIGN_VALIDATE.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT_VALIDATE.out.versions.first())

    SAMTOOLS_SORT_VALIDATE.out.bam.join(
        SAMTOOLS_INDEX_VALIDATE ( SAMTOOLS_SORT_VALIDATE.out.bam ).bai,
        by: [0], remainder: true)
        .set { ch_validate_bam_bai }


    ch_to_dupe_or_not = Channel.empty()
    // Check if REMOVE_DUPLICATES.out.bam exists.
    if (params.remove_duplicates) {
        REMOVE_DUPLICATES ( ch_validate_bam_bai )
        ch_versions = ch_versions.mix(REMOVE_DUPLICATES.out.versions.first())
        ch_to_dupe_or_not = REMOVE_DUPLICATES.out.bam_bai
    } else {
        ch_to_dupe_or_not = ch_validate_bam_bai
    }

    ch_to_dupe_or_not
        .join(EXTRACT_CHIMERIC_GENOMIC_TARGETS.out.gtf_extract, by: [0])
        .set { ch_validate_bam_bai_gtf }

    CHIMERIC_CONTIG_EVIDENCE_ANALYZER (
        ch_validate_bam_bai_gtf
    )
    ch_versions = ch_versions.mix(CHIMERIC_CONTIG_EVIDENCE_ANALYZER.out.versions.first())

    CHIMERIC_CONTIG_EVIDENCE_ANALYZER.out.evidence_bam
        .join(CHIMERIC_CONTIG_EVIDENCE_ANALYZER.out.evidence_bai, by: [0])
        .join(ABRIDGED_TSV.out.filtered_abridged, by: [0])
        .join(CHIMERIC_CONTIG_EVIDENCE_ANALYZER.out.evidence_counts, by: [0])
        .join(EXTRACT_CHIMERIC_GENOMIC_TARGETS.out.gtf_extract, by: [0])
        .join(EXTRACT_CHIMERIC_GENOMIC_TARGETS.out.fasta_extract, by: [0])
        .join(VIRUS_REPORT.out.genome_abundance_plot, by: [0])
        .join(VIRUS_REPORT.out.read_counts_image, by: [0])
        .join(VIRUS_REPORT.out.read_counts_log_image, by: [0])
        .set { ch_summary_report }

    SUMMARY_REPORT(
        ch_summary_report,
        params.gtf,
        ch_igvjs_VIF
    )
    ch_versions = ch_versions.mix(SUMMARY_REPORT.out.versions.first())

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
        //.mix(TRIMMOMATIC.out.mqc_log.collect{it[1]}.ifEmpty([]))
        .mix(PLUS.out.log_final.collect{it[1]}.ifEmpty([]))
        .mix(STAR_ALIGN_VALIDATE.out.log_final.collect{it[1]}.ifEmpty([]))


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
