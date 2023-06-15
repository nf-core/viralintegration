//
// Search for TE insertions
//

include { SAMTOOLS_VIEW } from '../../modules/nf-core/samtools/view/main.nf'
include { SAMTOOLS_FAIDX } from '../../modules/nf-core/samtools/faidx/main.nf'
include { BEDTOOLS_CLOSEST} from '../../modules/nf-core/bedtools/closest/main.nf'
include { BEDTOOLS_INTERSECT} from '../../modules/nf-core/bedtools/intersect/main.nf'
include { BEDTOOLS_MERGE as BEDTOOLS_MERGE_P
          BEDTOOLS_MERGE as BEDTOOLS_MERGE_M } from '../../modules/nf-core/bedtools/merge/main.nf'
include { BEDTOOLS_BAMTOBED} from '../../modules/nf-core/bedtools/bamtobed/main.nf'
include { PICARD_FILTERSAMREADS} from '../../modules/nf-core/picard/filtersamreads/main.nf'

include { DISCORDANT_READ_IDS } from '../../modules/local/discordant_read_ids.nf'
include { PRIMARY_READS_MERGED_P } from '../../modules/local/primary_reads_merged_p.nf'
include { PRIMARY_READS_MERGED_M } from '../../modules/local/primary_reads_merged_m.nf'


workflow TE_INSERTIONS {
    take:
    bam
    fasta
    te_gtf
    te_txt

    main:
    ch_versions = Channel.empty()

    // TODO transform gtf into individual TE GFF - probably a local module if we have to use the gff
    //  grep -P '[^(\w|\d|\-|\_|\#|\.)]'${line}'[^(\w|\d|\-|\_|\#|\.)]' $gtf > ${currdir}/${line}_TE.gff

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

    // TODO add conditional argument
    //if discordant readID file exists, then continue with remainder of TE analysis
    //then continue with following steps

    PICARD_FILTERSAMREADS (
        SAMTOOLS_VIEW.out.bam,
        DISCORDANT_READ_IDS.out.read_ids
    )
    ch_versions = ch_versions.mix(PICARD_FILTERSAMREADS.out.versions.first())

    BEDTOOLS_BAMTOBED (
        PICARD_FILTERSAMREADS.out.bam
    )
    ch_versions = ch_versions.mix(BEDTOOLS_BAMTOBED.out.versions.first())


    BEDTOOLS_MERGE_P (
        BEDTOOLS_BAMTOBED.out.bed
    )
    ch_versions = ch_versions.mix(BEDTOOLS_MERGE_P.out.versions.first())

    PRIMARY_READS_MERGED_P (
        BEDTOOLS_MERGE_P.out.bed
    )
    ch_versions = ch_versions.mix(PRIMARY_READS_MERGED_P.out.versions.first())

    BEDTOOLS_MERGE_M (
        BEDTOOLS_BAMTOBED.out.bed
    )
    ch_versions = ch_versions.mix(BEDTOOLS_MERGE_M.out.versions.first())

    PRIMARY_READS_MERGED_M (
        BEDTOOLS_MERGE_M.out.bed
    )
    ch_versions = ch_versions.mix(PRIMARY_READS_MERGED_M.out.versions.first())

    PRIMARY_READS_MERGE_P.out.plus_cluster
        .join(PRIMARY_READS_MERGE_M.out.minus_cluster, by: [0], remainder: true)
        .set { ch_plus_minus_bed }

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
