process ABRIDGED_TSV {
    tag "$meta.id"
    label 'process_low'

    // TODO Use python 3.6.9 and pigz in their own container
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "Conda environments cannot be used when using the PolyA-stripper script. Please use docker or singularity containers."
    }
    container "docker.io/trinityctat/ctat_vif"

    input:
    tuple val(meta), path(full_tsv)

    output:
    tuple val(meta), path("*.filtered.tsv")                , emit: filtered
    tuple val(meta), path("*.filtered.abridged.tsv")       , emit: filtered_abridged
    // TODO path "versions.yml"                  , emit: versions

    script: // This script is bundled with the pipeline, in nf-core/viralintegration/bin/
    def prefix = task.ext.prefix ?: "${meta.id}.vif.init"
    """
    #!/usr/bin/env python

    import pandas as pd
    min_reads = ${params.min_reads}

    # write abridged tsv
    df = pd.read_csv("$full_tsv", sep="\\t")
    df.drop('readnames', axis=1).to_csv("${prefix}.full.abridged.tsv", sep="\\t", index=False)

    #df = df[ (df.hits <= max_hits) & (df.total >= min_reads)]
    df = df[ df.total >= min_reads ]

    df.to_csv("${prefix}.filtered.tsv", sep="\\t", index=False)
    df.drop('readnames', axis=1).to_csv("${prefix}.filtered.abridged.tsv", sep="\\t", index=False)
    """
}
