//
// Align filtered reads against combined host and viral reference.
//

include { STAR_ALIGN as STAR_ALIGN_PLUS } from '../../modules/nf-core/star/align/main.nf'
include { STAR_GENOMEGENERATE } from '../../modules/nf-core/star/genomegenerate/main.nf'
include { SAMTOOLS_SORT } from '../../modules/nf-core/samtools/sort/main.nf'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main.nf'

workflow PLUS_TE {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    fasta // file: /path/to/fasta/
    gtf //   file: /path/to/gtf/

    main:
    ch_versions = Channel.empty()

    STAR_GENOMEGENERATE (
        fasta,
        gtf
    )
    ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions.first())
    ch_star_index = STAR_GENOMEGENERATE.out.index

    STAR_ALIGN_PLUS (
        reads,
        ch_star_index,
        gtf,
        false,
        "illumina",
        false
    )
    ch_versions = ch_versions.mix(STAR_ALIGN_PLUS.out.versions.first())

    SAMTOOLS_SORT (
        STAR_ALIGN_PLUS.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    SAMTOOLS_INDEX (
        SAMTOOLS_SORT.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    // TODO Concatenate files from HOST and PLUS
    // bam, bai, junction, sam_bam
    // Then join into channel

    SAMTOOLS_SORT.out.bam
        .join(SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
        .join(STAR_ALIGN_PLUS.out.junction)
        .set { ch_bam_bai_junction }


    emit:
    reads                  // channel: [ val(meta), [ reads ] ]

    index = ch_star_index

    // TODO Outputs will come from concatenation step
    // bam = STAR_ALIGN_PLUS.out.bam
    // junction = STAR_ALIGN_PLUS.out.junction
    log_final = STAR_ALIGN_PLUS.out.log_final
    log_out = STAR_ALIGN_PLUS.out.log_out
    log_progress = STAR_ALIGN_PLUS.out.log_progress

    //sam_bam = SAMTOOLS_SORT.out.bam

    //bai = SAMTOOLS_INDEX.out.bai

    bam_bai_junction = ch_bam_bai_junction

    versions = ch_versions // channel: [ versions.yml ]
}
