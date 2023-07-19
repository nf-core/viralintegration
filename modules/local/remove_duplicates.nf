process REMOVE_DUPLICATES {
    tag "$meta.id"
    label 'process_medium'

    // TODO Use python 3.6.9 and pigz in their own container
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "Conda environments cannot be used when using the PolyA-stripper script. Please use docker or singularity containers."
    }
    container "docker.io/trinityctat/ctat_vif"

    input:
    tuple val(meta), path(input_bam), path(input_bai)

    output:
    tuple val(meta), path("*.bam"), path("*.bai"), emit: bam_bai
    path "versions.yml"           , emit: versions

    script: // This script is bundled with the pipeline, in nf-core/viralintegration/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bam_mark_duplicates.py \
        -i ${input_bam} \
        -o ${prefix}.dedup.bam \
        -r

    samtools index ${prefix}.dedup.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
