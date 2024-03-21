process SUMMARY_REPORT {
    tag "$meta.id"
    label 'process_medium'

    // TODO Use python 3.6.9 and pigz in their own container
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "Conda environments cannot be used when using the PolyA-stripper script. Please use docker or singularity containers."
    }
    container "docker.io/trinityctat/ctat_vif"

    input:
    tuple val(meta),
        path(alignment_bam),
        path(alignment_bai),
        path(init_counts),
        path(vif_counts),
        path(chim_targets_gtf),
        path(chim_targets_fasta),
        path(genome_abundance_plot),
        path(read_counts_image),
        path(read_counts_log_image)
    path gtf
    path igvjs_VIF

    output:
    tuple val(meta), path ("*.html")                       , emit: html
    tuple val(meta), path ("*.prelim.refined.tsv")         , emit: prelim_refined_counts
    tuple val(meta), path ("*.refined.tsv")                , emit: refined_counts
    tuple val(meta), path ("*.refined.wRefGeneAnnots.tsv") , emit: refined_counts_w_genes
    tuple val(meta), path ("*.refined.distilled.tsv")      , emit: refined_distilled
    tuple val(meta), path ("*.genome_plot.png")            , emit: genome_abundance_plot
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when
    alignment_bam.size() > 0

    script: // This script is bundled with the pipeline, in nf-core/viralintegration/bin/
    // TODO Move to modules.config?
    def min_flank_frac_uniq = 0.0
    def max_coverage = 100
    def prefix = task.ext.prefix ?: "${meta.id}.vif"
    """
    refine_VIF_output.R \\
        --prelim_counts ${init_counts} \\
        --vif_counts ${vif_counts} \\
        --output ${prefix}.prelim.refined.tsv

    examine_flanking_uniq_kmer_composition.py \\
        --vif_tsv ${prefix}.prelim.refined.tsv \\
        --min_frac_uniq ${min_flank_frac_uniq} \\
        --output ${prefix}.refined.tsv

    distill_to_primary_target_list_via_brkpt_homologies.py \\
        --vif_tsv ${prefix}.refined.tsv \\
        > ${prefix}.refined.distilled.tsv

    make_VIF_genome_abundance_plot.R \\
        --vif_report ${prefix}.refined.tsv \\
        --title "Genome Wide Abundance" \\
        --output_png ${prefix}.genome_plot.png

    find_closest.py \\
        -i ${prefix}.refined.tsv \\
        -o ${prefix}.refined.wRefGeneAnnots.tsv \\
        --gtf ${gtf}

    if [[ -e ${prefix}.refined.wRefGeneAnnots.tsv ]]; then

        create_insertion_site_inspector_js.py \\
            --VIF_summary_tsv ${prefix}.refined.wRefGeneAnnots.tsv \\
            --json_outfile ${prefix}.igv.json

        # make bed for igvjs
        region_gtf_to_bed.py \\
            ${chim_targets_gtf} \\
            > ${prefix}.bed

        # prep for making the report
        # HACK Let's not hard code this in the future
        /usr/local/src/CTAT-VirusIntegrationFinder/util/bamsifter/bamsifter \\
            -c ${max_coverage} \\
            -o ${prefix}.reads.bam \\
            ${alignment_bam}

        # IGV reports expects to find, __PREFIX__.fa, __PREFIX__.bed, __PREFIX__.reads.bam
        ln -sf ${chim_targets_fasta} ${prefix}.fa

        make_VIF_igvjs_html.py \\
            --html_template $igvjs_VIF \\
            --fusions_json ${prefix}.igv.json \\
            --input_file_prefix ${prefix} \\
            --html_output ${prefix}.html

        # generate the final report
            add_to_html.py \\
            --html ${prefix}.html \\
            --out ${prefix}.html \\
            --image ${prefix}.genome_plot.png \\
            --image ${genome_abundance_plot} \\
            --image ${read_counts_image} \\
            --image ${read_counts_log_image}

    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
