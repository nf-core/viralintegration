//
// Refine chimeric read counts and create final summary report.
//

include { REMOVE_DUPLICATES } from '../../modules/local/remove_duplicates.nf'
include { CHIMERIC_CONTIG_EVIDENCE_ANALYZER } from '../../modules/local/chimeric_contig_evidence_analyzer.nf'
include { SUMMARY_REPORT } from '../../modules/local/summary_report.nf'

workflow CHIMERIC_READS {
    take:
    bam_bai
    gtf_extract
    fasta_extract
    filtered_abridged
    genome_abundance_plot
    read_counts_image
    read_counts_log_image
    gtf
    igvjs_VIF

    main:
    ch_versions = Channel.empty()

    ch_to_dupe_or_not = Channel.empty()
    // Check if REMOVE_DUPLICATES.out.bam exists.
    if (params.remove_duplicates) {
        REMOVE_DUPLICATES ( bam_bai )
        ch_versions = ch_versions.mix(REMOVE_DUPLICATES.out.versions.first())
        ch_to_dupe_or_not = REMOVE_DUPLICATES.out.bam_bai
    } else {
        ch_to_dupe_or_not = bam_bai
    }

    ch_to_dupe_or_not
        .join(gtf_extract, by: [0])
        .set { ch_validate_bam_bai_gtf }

    CHIMERIC_CONTIG_EVIDENCE_ANALYZER (
        ch_validate_bam_bai_gtf
    )
    ch_versions = ch_versions.mix(CHIMERIC_CONTIG_EVIDENCE_ANALYZER.out.versions.first())

    CHIMERIC_CONTIG_EVIDENCE_ANALYZER.out.evidence_bam
        .join(CHIMERIC_CONTIG_EVIDENCE_ANALYZER.out.evidence_bai, by: [0])
        .join(filtered_abridged, by: [0])
        .join(CHIMERIC_CONTIG_EVIDENCE_ANALYZER.out.evidence_counts, by: [0])
        .join(gtf_extract, by: [0])
        .join(fasta_extract, by: [0])
        .join(genome_abundance_plot, by: [0])
        .join(read_counts_image, by: [0])
        .join(read_counts_log_image, by: [0])
        .set { ch_summary_report }

    SUMMARY_REPORT(
        ch_summary_report,
        gtf,
        igvjs_VIF
    )
    ch_versions = ch_versions.mix(SUMMARY_REPORT.out.versions.first())


    emit:

    evidence_counts = CHIMERIC_CONTIG_EVIDENCE_ANALYZER.out.evidence_counts

    html = SUMMARY_REPORT.out.html
    refined_counts = SUMMARY_REPORT.out.refined_counts
    refined_counts_w_genes = SUMMARY_REPORT.out.refined_counts_w_genes
    genome_abundance_plot = SUMMARY_REPORT.out.genome_abundance_plot

    versions = ch_versions // channel: [ versions.yml ]
}
