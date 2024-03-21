process VIRUS_REPORT {
    tag "$meta.id"
    label 'process_medium'

    // TODO Use python 3.6.9 and pigz in their own container
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "Conda environments cannot be used when using the PolyA-stripper script. Please use docker or singularity containers."
    }
    container "docker.io/trinityctat/ctat_vif"

    input:
    tuple val(meta), path(bam), path(bai), path(insertion_site_candidates)
    path viral_fasta
    path igvjs_VIF

    output:
    tuple val(meta), path ("*.igvjs.html")                    , optional: true, emit: html
    tuple val(meta), path ("*.init.genome_plot.png")          , emit: genome_abundance_plot
    tuple val(meta), path ("*.igvjs.bam")                     , emit: virus_alignments_bam
    tuple val(meta), path ("*.igvjs.bam.bai")                 , emit: virus_alignments_bai
    tuple val(meta), path ("*.virus_read_counts_summary.tsv") , optional: true, emit: read_counts_summary
    tuple val(meta), path ("*.virus_read_counts.png")         , emit: read_counts_image
    tuple val(meta), path ("*.virus_read_counts_log.png")     , emit: read_counts_log_image
    tuple val(meta), path ("*.virus_coverage_*.png")          , optional: true, emit: virus_images
    path "versions.yml"                                       , emit: versions

    script: // This script is bundled with the pipeline, in nf-core/viralintegration/bin/
    // TODO Move to modules.config?
    def prefix = task.ext.prefix ?: "${meta.id}.VirusDetect"
    def remove_duplicates = '--remove_duplicates'
    // TODO remove duplicates option
    def num_top_viruses = '--num_top_viruses 20'
    """
    make_VIF_genome_abundance_plot.R \\
        --vif_report ${insertion_site_candidates} \\
        --title "Preliminary Genome Wide Abundance" \\
        --output_png ${prefix}.init.genome_plot.png

    # restrict bam to only viruses of interest
    bam=${bam}
    samtools faidx ${viral_fasta}
    awk '{printf("%s\t0\t%s\\n",\$1,\$2);}' ${viral_fasta}.fai  > viruses.bed
    samtools view -b -L viruses.bed \$bam -o ${prefix}.igvjs.bam
    samtools sort -o ${prefix}.sorted.igvjs.bam -T $prefix ${prefix}.igvjs.bam
    bam="${prefix}.sorted.igvjs.bam"
    samtools index \$bam

    # clean up the bam, restrict to proper pairs and non-supplemental alignments
    restrict_bam_to_proper_aligns.py \$bam \${bam}.clean.bam
    mv \${bam}.clean.bam \${bam}
    mv \${bam}.clean.bam.bai \${bam}.bai


    if [ "${remove_duplicates}" == "true" ]; then
        bam_mark_duplicates.py -i ${bam} -o dups.removed.bam -r
        mv dups.removed.bam  ${bam}
        samtools index ${bam}
    fi

    # generates read_counts_summary and images
    plot_top_virus_coverage.R \\
        --vif_report ${insertion_site_candidates}  \\
        --virus_fai ${viral_fasta}.fai \\
        --bam \${bam} \\
        --output_prefix ${prefix}

    if [[ -s "${prefix}.virus_read_counts_summary.tsv" ]] ; then
        # make bed for igvjs
        create_igvjs_virus_bed.py \\
            --summary ${prefix}.virus_read_counts_summary.tsv \\
            --output_prefix ${prefix} \\
            ${num_top_viruses}

        create_insertion_site_inspector_js.py \\
            --VIF_summary_tsv ${prefix}.igvjs.table.tsv \\
            --json_outfile ${prefix}.igvjs.json

        # prep for making the report
        # HACK Let's not hard code this in the future
        /usr/local/src/CTAT-VirusIntegrationFinder/util/bamsifter/bamsifter \\
            -c 100 \\
            -o ${prefix}.igvjs.reads.bam \\
            \${bam}

        # IGV reports expects to find, __PREFIX__.fa, __PREFIX__.bed, __PREFIX__.reads.bam
        #ln -sf ${viral_fasta} ${prefix}.virus.fa
        create_igvjs_virus_fa.py \\
            ${prefix}.igvjs.bed \\
            ${viral_fasta}  \\
            ${prefix}.igvjs.fa

        # generate the html
        make_VIF_igvjs_html.py \\
                --html_template $igvjs_VIF \\
                --fusions_json ${prefix}.igvjs.json \\
                --input_file_prefix ${prefix}.igvjs \\
                --html_output ${prefix}.igvjs.html
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
