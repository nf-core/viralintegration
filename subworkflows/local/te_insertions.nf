//
// Search for TE insertions
//

include { SAMTOOLS_VIEW } from '../../modules/nf-core/samtools/view/main.nf'
include { SAMTOOLS_FAIDX } from '../../modules/nf-core/samtools/faidx/main.nf'
include { BEDTOOLS_CLOSEST} from '../../modules/nf-core/bedtools/closest/main.nf'
include { BEDTOOLS_INTERSECT} from '../../modules/nf-core/bedtools/intersect/main.nf'
include { BEDTOOLS_MERGE} from '../../modules/nf-core/bedtools/merge/main.nf'
include { BEDTOOLS_BAMTOBED} from '../../modules/nf-core/bedtools/bamtobed/main.nf'
include { PICARD_FILTERSAMREADS} from '../../modules/nf-core/picard/filtersamreads/main.nf'
include { DISCORDANT_READ_IDS } from '../../modules/local/discordant_read_ids.nf'


workflow TE_INSERTIONS {
    take:
    bam
    fasta
    te_gtf
    te_txt

    main:
    ch_versions = Channel.empty()

    bam
        .join(te_gtf, by: [0], remainder: true)
        .set { ch_bed_intersect }

    BEDTOOLS_INTERSECT (
        ch_bed_intersect,
        "filter.bed" // TODO Do we need a chromosome size file? - this may be the te_txt list with TEs of interest
    )
    ch_versions = ch_versions.mix(BEDTOOLS_INTERSECT.out.versions.first())

    SAMTOOLS_VIEW (
        BEDTOOLS_INTERSECT.out.bam,
        fasta
    )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions.first())

    DISCORDANT_READ_IDS (
        SAMTOOLS_VIEW.out.bam
    )
    ch_versions = ch_versions.mix(DISCORDANT_READ_IDS.out.versions.first())

    //if discordant readID file exists, then continue with remainder of TE analysis
    //then continue with following steps

    PICARD_FILTERSAMREADS (
        SAMTOOLS_VIEW.out.bam,
        DISCORDANT_READ_IDS.out.read_ids
    )
    ch_versions = ch_versions.mix(PICARD_FILTERSAMREADS.out.versions.first())


    // awk '{if ($4 > 3) print $0}' > ${currdir}/${line}_plusCluster.bed


    BEDTOOLS_BAMTOBED (
        PICARD_FILTERSAMREADS.out.bam
    )
    ch_versions = ch_versions.mix(BEDTOOLS_BAMTOBED.out.versions.first())


    BEDTOOLS_MERGE (
        BEDTOOLS_BAMTOBED.out.bed
    )
    ch_versions = ch_versions.mix(BEDTOOLS_MERGE.out.versions.first())

    // TODO two steps with plus and minus strand - do we need 2 bed files?
    // awk '{if ($4 > 3) print $0}' > ${currdir}/${line}_plusCluster.bed
    // awk '{if ($4 > 3) print $0}' > ${currdir}/${line}_minusCluster.bed

 // TODO make a channel for tupule for plus and minus
    // plusCluster.bed
    //     .join(minusCluster.bed, by: [0], remainder: true)
    //     .set { ch_plus_minus_bed }

    SAMTOOLS_FAIDX(
        fasta,
        fai // TODO do we need this?
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions.first())


    BEDTOOLS_CLOSEST (
        ch_plus_minus_bed,
        SAMTOOLS_FAIDX.out.fai
    )
    ch_versions = ch_versions.mix(BEDTOOLS_CLOSEST.out.versions.first())



    // awk stuff

    emit:

    //id emitting = calling the file
    bam_bai = ch_validate_bam_bai
    log_final = STAR_ALIGN_VALIDATE.out.log_final

    versions = ch_versions // channel: [ versions.yml ]
}
