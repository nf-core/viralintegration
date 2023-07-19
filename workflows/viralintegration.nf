/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowViralintegration.initialise(params, log)

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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { TRIMMOMATIC                 } from '../modules/nf-core/trimmomatic/main'
include { STAR_GENOMEGENERATE as STAR_GENOMEGENERATE_HOST
          STAR_GENOMEGENERATE as STAR_GENOMEGENERATE_PLUS } from '../modules/nf-core/star/genomegenerate/main'
include { STAR_ALIGN as STAR_ALIGN_HOST
          STAR_ALIGN as STAR_ALIGN_PLUS } from '../modules/nf-core/star/align/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_PLUS
          SAMTOOLS_SORT as SAMTOOLS_SORT_VALIDATE
          SAMTOOLS_SORT as SAMTOOLS_SORT_DUPLICATES } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_PLUS
          SAMTOOLS_INDEX as SAMTOOLS_INDEX_VALIDATE
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
        file(params.input)
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    STAR_GENOMEGENERATE_HOST (
        params.fasta,
        params.gtf
    )
    ch_versions = ch_versions.mix(STAR_GENOMEGENERATE_HOST.out.versions)

    // TODO Use igenomes
    STAR_ALIGN_HOST (
        INPUT_CHECK.out.reads,
        STAR_GENOMEGENERATE_HOST.out.index,
        params.gtf,
        false,
        "illumina",
        false
    )
    ch_versions = ch_versions.mix(STAR_ALIGN_HOST.out.versions.first())

    TRIMMOMATIC (
        STAR_ALIGN_HOST.out.fastq
    )
    ch_versions = ch_versions.mix(TRIMMOMATIC.out.versions.first())

    POLYA_STRIPPER (
        TRIMMOMATIC.out.trimmed_reads
    )
    ch_versions = ch_versions.mix(POLYA_STRIPPER.out.versions.first())

    CAT_FASTA (
        params.fasta,
        params.viral_fasta
    )
    ch_versions = ch_versions.mix(CAT_FASTA.out.versions)

    STAR_GENOMEGENERATE_PLUS (
        CAT_FASTA.out.plus_fasta,
        params.gtf
    )
    ch_versions = ch_versions.mix(STAR_GENOMEGENERATE_PLUS.out.versions)

    STAR_ALIGN_PLUS (
        POLYA_STRIPPER.out.polya_trimmed,
        STAR_GENOMEGENERATE_PLUS.out.index,
        params.gtf,
        false,
        "illumina",
        false
    )
    ch_versions = ch_versions.mix(STAR_ALIGN_PLUS.out.versions.first())

    SAMTOOLS_SORT_PLUS (
        STAR_ALIGN_PLUS.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT_PLUS.out.versions.first())

    SAMTOOLS_INDEX_PLUS (
        SAMTOOLS_SORT_PLUS.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_PLUS.out.versions.first())

    SAMTOOLS_SORT_PLUS.out.bam
        .join(SAMTOOLS_INDEX_PLUS.out.bai, by: [0], remainder: true)
        .join(STAR_ALIGN_PLUS.out.junction)
        .set { ch_bam_bai_junction }

    INSERTION_SITE_CANDIDATES (
        ch_bam_bai_junction,
        params.fasta,
        params.viral_fasta
    )
    ch_versions = ch_versions.mix(INSERTION_SITE_CANDIDATES.out.versions.first())

    ABRIDGED_TSV (
        INSERTION_SITE_CANDIDATES.out.full
    )
    // TODO ch_versions = ch_versions.mix(ABRIDGED_TSV.out.versions.first())

    SAMTOOLS_SORT_PLUS.out.bam
        .join(SAMTOOLS_INDEX_PLUS.out.bai, by: [0], remainder: true)
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

    STAR_ALIGN_HOST.out.fastq
        .join(EXTRACT_CHIMERIC_GENOMIC_TARGETS.out.fasta_extract)
        .set { ch_unaligned_fastq_fasta }
    STAR_ALIGN_VALIDATE (
        ch_unaligned_fastq_fasta,
        STAR_GENOMEGENERATE_PLUS.out.index,
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

    methods_description    = WorkflowViralintegration.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
        .mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        .mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
        .mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
        .mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
        .mix(STAR_ALIGN_HOST.out.log_final.collect{it[1]}.ifEmpty([]))
        .mix(TRIMMOMATIC.out.mqc_log.collect{it[1]}.ifEmpty([]))
        .mix(STAR_ALIGN_PLUS.out.log_final.collect{it[1]}.ifEmpty([]))
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
