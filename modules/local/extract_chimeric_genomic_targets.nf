process EXTRACT_CHIMERIC_GENOMIC_TARGETS {
    tag "$insertion_site_candidates_abridged"

    // TODO Use python 3.6.9 and pigz in their own container
    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the PolyA-stripper script. Please use docker or singularity containers."
    }
    container "trinityctat/ctat_vif"

    input:
    path fasta
    path viral_fasta
    path insertion_site_candidates_abridged

    output:
    path ("*.fasta")                           , emit: fasta_extract
    path ("*.gtf")                             , emit: gtf_extract
    path "versions.yml"                        , emit: versions

    script: // This script is bundled with the pipeline, in nf-core/viralintegration/bin/
    // TODO Move to modules.config?
    def prefix = task.ext.prefix ?: "${insertion_site_candidates_abridged}.vif.extract"
    def pad_region_length = '--pad_region_length 1000'
    """
    extract_chimeric_genomic_targets.py \\
        --fasta ${fasta} \\
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