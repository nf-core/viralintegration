process CAT_JUNCTION {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    tuple val(meta), path (host_chim_junction)
    tuple val(meta), path (plus_chim_junction)

    output:
    tuple val(meta), path ("*_chim.junction") , emit: chim_junction
    path "versions.yml"                       , emit: versions

    script:
    """
    sed -i '1d' $plus_chim_junction

    cat $host_chim_junction $plus_chim_junction > ${meta.id}_chim.junction

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
    END_VERSIONS
    """
}
