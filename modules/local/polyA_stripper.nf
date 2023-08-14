process POLYA_STRIPPER {
    tag "$meta.id"
    label 'process_medium'

    // TODO Use python 3.6.9 and pigz in their own container
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "Conda environments cannot be used when using the PolyA-stripper script. Please use docker or singularity containers."
    }
    container "docker.io/trinityctat/ctat_vif"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.polyA-trimmed.fastq.gz'), emit: polya_trimmed
    path "versions.yml"                              , emit: versions

    script: // This script is bundled with the pipeline, in nf-core/viralintegration/bin/
    """
    fastq_polyA_stripper.py \\
        --out_prefix ${meta.id} \\
        --left_fq ${reads[0]} \\
    ### --right_fq ${reads[1]} # HACK for single end reads, TODO make this input optional

    pigz *.polyA-trimmed.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
