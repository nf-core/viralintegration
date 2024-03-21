process EXTRACT_CHIMERIC_GENOMIC_TARGETS {
    tag "$meta.id"
    label 'process_low'

    // TODO Use python 3.6.9 and pigz in their own container
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "Conda environments cannot be used when using the PolyA-stripper script. Please use docker or singularity containers."
    }
    container "docker.io/trinityctat/ctat_vif"

    input:
    tuple val(meta), path(insertion_site_candidates_abridged)
    path(ref_genome_fasta)
    path(viral_fasta)

    output:
    tuple val(meta), path ("*.fasta")          , emit: fasta_extract
    tuple val(meta), path ("*.gtf")            , emit: gtf_extract
    path "versions.yml"                        , emit: versions

    script: // This script is bundled with the pipeline, in nf-core/viralintegration/bin/
    // TODO Move to modules.config?
    def prefix = task.ext.prefix ?: "${meta.id}.vif.extract"
    def pad_region_length = '--pad_region_length 1000'
    """
    extract_chimeric_genomic_targets.py \\
        --fasta ${ref_genome_fasta} \\
        --patch_db_fasta ${viral_fasta} \\
        --output_prefix ${prefix} \\
        --chim_events ${insertion_site_candidates_abridged} \\
        ${pad_region_length}


    if [ -s ${prefix}.fasta ]; then echo "true" > has_results; else echo "false" > has_results; fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
