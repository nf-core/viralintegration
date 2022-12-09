process CAT_FASTA {
    tag "$ref_genome_fasta"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    path ref_genome_fasta
    path virus_db_fasta

    output:
    path '*_plus_viraldb.fasta', emit: plus_fasta
    path "versions.yml"        , emit: versions

    script:
    """
    cat $ref_genome_fasta $virus_db_fasta > ${ref_genome_fasta.simpleName}_plus_viraldb.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
    END_VERSIONS
    """
}
