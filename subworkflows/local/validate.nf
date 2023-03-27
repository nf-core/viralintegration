//
// Validate potential chimeric reads.
//

include { STAR_ALIGN_VALIDATE } from '../../modules/local/star_align_validate.nf'
include { SAMTOOLS_SORT} from '../../modules/nf-core/samtools/sort/main.nf'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main.nf'

workflow VALIDATE {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    fasta_extract
    index

    main:
    ch_versions = Channel.empty()

        reads
        .join(fasta_extract)
        .set { ch_unaligned_fastq_fasta }

    STAR_ALIGN_VALIDATE (
        ch_unaligned_fastq_fasta,
        index,
        "illumina",
        false
    )
    ch_versions = ch_versions.mix(STAR_ALIGN_VALIDATE.out.versions.first())

    SAMTOOLS_SORT (
        STAR_ALIGN_VALIDATE.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    SAMTOOLS_SORT.out.bam.join(
        SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam ).bai,
        by: [0], remainder: true)
        .set { ch_validate_bam_bai }


    emit:

    bam_bai = ch_validate_bam_bai

    versions = ch_versions // channel: [ versions.yml ]
}
