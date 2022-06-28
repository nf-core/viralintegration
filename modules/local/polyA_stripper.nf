process POLYA_STRIPPER {
    tag "$meta.id"

    conda (params.enable_conda ? "conda-forge::python=3.6.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        // HACK The script breaks under 3.8.3
        'python:3.6.9' :
        'python:3.6.9' }"

    input:
    tuple val(meta), path(reads)

    output:
    path '*.polyA-trimmed.fastq', emit: polya_trimmed
    path "versions.yml"         , emit: versions

    script: // This script is bundled with the pipeline, in nf-core/viralintegration/bin/
    """
    fastq_polyA_stripper.py \\
        --out_prefix ${meta.id} \\
        --left_fq ${reads[0]} \\
        --right_fq ${reads[1]}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
