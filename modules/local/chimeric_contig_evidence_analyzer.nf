process CHIMERIC_CONTIG_EVIDENCE_ANALYZER {
    tag "$meta.id"
    label 'process_low'

    // TODO Use python 3.6.9 and pigz in their own container
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "Conda environments cannot be used when using the PolyA-stripper script. Please use docker or singularity containers."
    }
    container "docker.io/trinityctat/ctat_vif"

    input:
    tuple val(meta), path(bam), path(bai), path(gtf)

    output:
    tuple val(meta), path ("*.evidence_counts.tsv")   , emit: evidence_counts
    tuple val(meta), path ("*.evidence.bam")          , emit: evidence_bam
    tuple val(meta), path ("*.evidence.bam.bai")      , emit: evidence_bai
    tuple val(meta), path ("insertions_validated.txt"), emit: txt
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when
    gtf.size() > 0

    script: // This script is bundled with the pipeline, in nf-core/viralintegration/bin/
    // TODO Move to modules.config?
    def prefix = task.ext.prefix ?: "${meta.id}.vif"
    def insertions_file = "${prefix}.evidence_counts.tsv"
    def insertions_validated_file= "insertions_validated.txt"
    def num_insertions= "(cat ${insertions_file} | wc -l)"

    // TODO Fix insertions_validated output
    """
    chimeric_contig_evidence_analyzer.py \
        --patch_db_bam ${bam} \
        --patch_db_gtf ${gtf} \
        --output_prefix ${prefix}

    samtools index ${prefix}.evidence.bam

    if (( ${num_insertions} > 1 )); then
        echo "true" > ${insertions_validated_file}
    else
        echo "false" > ${insertions_validated_file}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
