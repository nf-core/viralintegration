process CAT_JUNCTION {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    path host_chim_junction
    path plus_chim_junction

    output:
    path '*_plus_chim.junction', emit: plus_chim_junction
    path "versions.yml"        , emit: versions

    script:
    """
    sed -i '1d' $plus_chim_junction

    cat $host_chim_junction $plus_chim_junction > ${meta.id}_plus_chim.junction

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
    END_VERSIONS
    """
}
