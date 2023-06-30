process CAT_JUNCTION {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    tuple val(meta), path (host_chim_junction), path (plus_chim_junction)

    output:
    tuple val(meta), path '*_chim.junction', emit: plus_chim_junction
    path "versions.yml"        , emit: versions

    script:
    def junction_file = "${meta.id}.*.Chimeric.out.junction"
    def line_num = "(cat ${junction_file} | wc -l)"
    """
    if ${line_num} > 1; then

        sed -i '1d' ${meta.id}plus_chim_junction

        cat ${meta.id}host_chim_junction ${meta.id}plus_chim_junction > ${meta.id}_chim.junction

    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
    END_VERSIONS
    """
}
