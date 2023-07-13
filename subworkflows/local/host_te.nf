//
// Align input reads against host genome.
//

include { TRIMMOMATIC } from "../../modules/nf-core/trimmomatic/main.nf"
include { POLYA_STRIPPER } from '../../modules/local/polyA_stripper.nf'

workflow HOST_TE {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()

    TRIMMOMATIC (
        reads
    )
    ch_versions = ch_versions.mix(TRIMMOMATIC.out.versions.first())

    POLYA_STRIPPER (
        TRIMMOMATIC.out.trimmed_reads
    )
    ch_versions = ch_versions.mix(POLYA_STRIPPER.out.versions.first())


    emit:
    reads                  // channel: [ val(meta), [ reads ] ]

    trimmed_reads = TRIMMOMATIC.out.trimmed_reads
    //mqc_log = TRIMMOMATIC.out.mqc_log

    polya_trimmed = POLYA_STRIPPER.out.polya_trimmed

    versions = ch_versions // channel: [ versions.yml ]
}
