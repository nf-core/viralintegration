//
// Align input reads against host genome.
//

include { STAR_ALIGN } from '../../modules/nf-core/star/align/main.nf'
include { STAR_GENOMEGENERATE_GENERATE } from '../../modules/nf-core/star/genomegenerate/main.nf'

workflow HOST_STAR {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    fasta // file: /path/to/fasta/
    gtf //   file: /path/to/gtf/

    main:
    STAR_GENOMEGENEREATE (
        fasta,
        gtf
    )
    ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions.first())
    ch_star_index = STAR_GENOMEGENERATE.out.index

    STAR_ALIGN (
        reads,
        ch_star_index,
        gtf,
        false,
        "illumina",
        false
    )
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())


    emit:
    reads                  // channel: [ val(meta), [ reads ] ]
    index = ch_star_index
    bam = STAR_ALIGN.out.bam
    log_final = STAR_ALIGN.out.log_final
    log_out = STAR_ALIGN.out.log_out
    log_progress = STAR_ALIGN.out.log_progress
    fastq = STAR_ALIGN.out.fastq
    versions = ch_versions // channel: [ versions.yml ]
}
