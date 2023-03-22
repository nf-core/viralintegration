//
// Map potential virus insertion sites and create virus infection evidence viewer.
//

include { INSERTION_SITE_CANDIDATES } from '../../modules/local/insertion_site_candidates.nf'
include { ABRIDGED_TSV } from '../../modules/local/abridged_tsv.nf'
include { VIRUS_REPORT } from '../../modules/local/virus_report.nf'
include { EXTRACT_CHIMERIC_GENOMIC_TARGETS } from '../../modules/local/extract_chimeric_genomic_targets.nf'

workflow INSERTION_SITES {
    take:
    junction // file: /path/to/junction/
    fasta // file: /path/to/fasta/
    viral_fasta // file: /path/to/viral_fasta/
    bam //   file: /path/to/bam/
    bai //   file: /path/to/bai/
    igvjs_VIF


    main:
    ch_versions = Channel.empty()

    INSERTION_SITE_CANDIDATES (
        junction,
        fasta,
        viral_fasta
    )
    ch_versions = ch_versions.mix(INSERTION_SITE_CANDIDATES.out.versions.first())

    ABRIDGED_TSV (
        INSERTION_SITE_CANDIDATES.out.full
    )
    // TODO ch_versions = ch_versions.mix(ABRIDGED_TSV.out.versions.first())

    bam
        .join(bai, by: [0], remainder: true)
        .join(ABRIDGED_TSV.out.filtered_abridged)
        .set { ch_bam_bai_filtered }

    VIRUS_REPORT (
        ch_bam_bai_filtered,
        viral_fasta,
        igvjs_VIF
    )
    ch_versions = ch_versions.mix(VIRUS_REPORT.out.versions.first())

    // TODO Handle insertion_site_candidates
    // File insertion_site_candidates_use = select_first([insertion_site_candidates, InsertionSiteCandidates.filtered_abridged])

    EXTRACT_CHIMERIC_GENOMIC_TARGETS (
        ABRIDGED_TSV.out.filtered_abridged,
        fasta,
        viral_fasta
    )
    ch_versions = ch_versions.mix(EXTRACT_CHIMERIC_GENOMIC_TARGETS.out.versions.first())

    emit:
    full = INSERTION_SITE_CANDIDATES.out.full
    genome_chimeric_evidence_reads_bam = INSERTION_SITE_CANDIDATES.out.genome_chimeric_evidence_reads_bam
    genome_chimeric_evidence_reads_bai = INSERTION_SITE_CANDIDATES.out.genome_chimeric_evidence_reads_bai

    filtered_abridged = ABRIDGED_TSV.out.filtered_abridged

    html = VIRUS_REPORT.out.html
    genome_abundance_plot = VIRUS_REPORT.out.genome_abundance_plot
    virus_alignments_bam = VIRUS_REPORT.out.virus_alignments_bam
    virus_alignments_bai = VIRUS_REPORT.out.virus_alignments_bai
    read_counts_summary = VIRUS_REPORT.out.read_counts_summary
    read_counts_image = VIRUS_REPORT.out.read_counts_image
    read_counts_log_image = VIRUS_REPORT.out.read_counts_log_image

    fasta_extract = EXTRACT_CHIMERIC_GENOMIC_TARGETS.out.fasta_extract
    gtf_extract = EXTRACT_CHIMERIC_GENOMIC_TARGETS.out.gtf_extract

    versions = ch_versions // channel: [ versions.yml ]
}
