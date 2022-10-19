/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowViralintegration.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()
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
include { FASTQC                      } from '../modules/nf-core/modules/fastqc/main'
include { TRIMMOMATIC                 } from '../modules/nf-core/modules/trimmomatic/main'
include { SAMTOOLS_FAIDX              } from '../modules/nf-core/modules/samtools/faidx/main'
include { STAR_GENOMEGENERATE as STAR_GENOMEGENERATE_HOST
          STAR_GENOMEGENERATE as STAR_GENOMEGENERATE_PLUS } from '../modules/nf-core/modules/star/genomegenerate/main'
include { STAR_ALIGN as STAR_ALIGN_HOST
          STAR_ALIGN as STAR_ALIGN_PLUS } from '../modules/nf-core/modules/star/align/main'
include { SAMTOOLS_SORT               } from '../modules/nf-core/modules/samtools/sort/main'
include { SAMTOOLS_INDEX              } from '../modules/nf-core/modules/samtools/index/main'
include { MULTIQC                     } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

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

    STAR_GENOMEGENERATE_HOST (
        params.fasta,
        params.gtf
    )

    // TODO Use igenomes
    STAR_ALIGN_HOST (
        INPUT_CHECK.out.reads,
        STAR_GENOMEGENERATE_HOST.out.index,
        params.gtf,
        false,
        "illumina",
        false
    )

    TRIMMOMATIC (
        STAR_ALIGN_HOST.out.fastq
    )

    POLYA_STRIPPER (
        TRIMMOMATIC.out.trimmed_reads
    )

    viral_fasta = [ [ id:'viral_fasta', single_end:false ], // meta map
                    file(params.viral_fasta, checkIfExists: true) ]

    CAT_FASTA (
        params.fasta,
        params.viral_fasta
    )

    SAMTOOLS_FAIDX (
        viral_fasta
    )

    STAR_GENOMEGENERATE_PLUS (
        CAT_FASTA.out.plus_fasta,
        params.gtf
    )

    STAR_ALIGN_PLUS (
        POLYA_STRIPPER.out.polya_trimmed,
        STAR_GENOMEGENERATE_PLUS.out.index,
        params.gtf,
        false,
        "illumina",
        false
    )

    SAMTOOLS_SORT (
        STAR_ALIGN_PLUS.out.bam
    )

    SAMTOOLS_INDEX (
        SAMTOOLS_SORT.out.bam
    )

    SAMTOOLS_SORT.out.bam
        .join(SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
        .join(STAR_ALIGN_PLUS.out.junction)
        .set { ch_bam_bai_junction }

    INSERTION_SITE_CANDIDATES (
        ch_bam_bai_junction,
        params.fasta,
        params.viral_fasta
    )

    ABRIDGED_TSV (
        INSERTION_SITE_CANDIDATES.out.full
    )

    SAMTOOLS_SORT.out.bam
        .join(SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
        .join(ABRIDGED_TSV.out.filtered_abridged)
        .set { ch_bam_bai_filtered }

    VIRUS_REPORT (
        ch_bam_bai_filtered,
        params.viral_fasta,
        ch_igvjs_VIF
    )

    // Join and pass fai and fasta pairs together //
    ch_fasta = Channel.of([[ id:'fasta', single_end:false ], // meta map
                           file(params.fasta, checkIfExists: true)])
        .view()

    ch_viral_fasta = Channel.of([[ id:'viral_fasta', single_end:false ], // meta map
                                 file(params.viral_fasta, checkIfExists: true)])
        .view()

    STAR_GENOMEGENERATE_HOST.out.index
        .join(ch_fasta, by: [0], remainder: true)
        .set { ch_ref_fa_fai }

    SAMTOOLS_FAIDX.out.fai
        .join(ch_viral_fasta, by: [0], remainder: true)
        .set { ch_viral_fa_fai }


    EXTRACT_CHIMERIC_GENOMIC_TARGETS (
        ABRIDGED_TSV.out.filtered_abridged,
        ch_ref_fa_fai,
        ch_viral_fa_fai
    )

    STAR_ALIGN_HOST.out.fastq
        .join(EXTRACT_CHIMERIC_GENOMIC_TARGETS.out.fasta_extract)
        .set { ch_unaligned_fastq_fasta }
    STAR_ALIGN_VALIDATE (
        ch_unaligned_fastq_fasta,
        STAR_GENOMEGENERATE_PLUS.out.index,
        "illumina",
        false
    )

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowViralintegration.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(TRIMMOMATIC.out.log.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
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
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
