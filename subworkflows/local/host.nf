//
// Align input reads against host genome.
//

include { STAR_ALIGN } from '../../modules/nf-core/star/align/main.nf'
include { STAR_GENOMEGENERATE } from '../../modules/nf-core/star/genomegenerate/main.nf'
include { TRIMMOMATIC } from "../../modules/nf-core/trimmomatic/main.nf"
include { POLYA_STRIPPER } from '../../modules/local/polyA_stripper.nf'

workflow HOST {
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

    STAR_ALIGN (
        reads,
        ch_star_index,
        gtf,
        false,
        "illumina",
        false
    )
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())

    TRIMMOMATIC (
        STAR_ALIGN.out.fastq
    )
    ch_versions = ch_versions.mix(TRIMMOMATIC.out.versions.first())

    POLYA_STRIPPER (
        TRIMMOMATIC.out.trimmed_reads
    )
    ch_versions = ch_versions.mix(POLYA_STRIPPER.out.versions.first())


    emit:
    reads                  // channel: [ val(meta), [ reads ] ]

    index = ch_star_index

    bam = STAR_ALIGN.out.bam
    log_final = STAR_ALIGN.out.log_final
    log_out = STAR_ALIGN.out.log_out
    log_progress = STAR_ALIGN.out.log_progress
    fastq = STAR_ALIGN.out.fastq

    trimmed_reads = TRIMMOMATIC.out.trimmed_reads

    polya_trimmed = POLYA_STRIPPER.out.polya_trimmed

    versions = ch_versions // channel: [ versions.yml ]
}
