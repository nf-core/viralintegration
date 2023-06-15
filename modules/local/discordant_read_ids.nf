process DISCORDANT_READ_IDS {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*_ReadID.txt"), emit: read_id
    path "versions.yml"                  , emit: versions

    script:
    def fis = 400
    """
    awk -v Ins=`expr ${fis} \* 10` '{if ((\$7 != "=") || (\$9 > Ins) || (\$9 < -Ins)) print \$1}' > \${line}_ReadID.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
    END_VERSIONS
    """
}
