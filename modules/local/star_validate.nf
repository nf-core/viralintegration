process STAR_VALIDATE {
    tag "$sample_id"

    // TODO Use python 3.6.9 and pigz in their own container
    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the PolyA-stripper script. Please use docker or singularity containers."
    }
    container "trinityctat/ctat_vif"

    input:
    path fastq1
    path fastq2
    path star_reference

    output:
    path ("*.Aligned.sortedByCoord.out.bam")     , emit: bam
    path ("*.Aligned.sortedByCoord.out.bam.bai") , emit: bai
    path ("*.Log.final.out")                     , emit: output_log_final
    path ("*.SJ.out.tab")                        , emit: output_SJ
    path ("*.Chimeric.out.junction")             , emit: chimeric_junction
    path ("versions.yml")                        , emit: versions

    script: // This script is bundled with the pipeline, in nf-core/viralintegration/bin/
    // TODO Move to modules.config?
    def prefix = task.ext.prefix ?: "${sample_id}.validate_inserts"
    def max_mate_dist = '100000'
    """
    cpu=~{cpu}
        genomeDir="~{star_reference_dirpath}"
        if [[ "${genomeDir}" == "" ]]; then
            genomeDir="~{star_reference}"
        fi

        fastqs="${fastq1} ${fastq2}"
        if [[ "${fastq2}" == "" ]] || [[ ! -s "${fastq2}" ]] ; then
        fastqs="${fastq1}"
        fi



        readFilesCommand=""
        if [[ "${fastq1}" == *.gz ]] ; then
            readFilesCommand="--readFilesCommand \"gunzip -c\""
        fi
        if [ "~{autodetect_cpu}" == "true" ]; then
            cpu=~{nproc}
        fi


        if [ -f "${genomeDir}" ] ; then
            mkdir genome_dir
            compress=""
            if [[ $genomeDir == *.tar.gz ]] ; then
                compress="-I pigz"
            elif [[ $genomeDir == *.tar.bz2 ]] ; then
                compress="-I pbzip2"
            fi
            tar $compress -xf ~{star_reference} -C genome_dir --strip-components 1
            genomeDir="genome_dir"
        fi

        # special case for tar of fastq files
        if [[ "${fastq1}" == *.tar.gz ]] ; then
            mkdir fastq
            tar -I pigz -xvf ${fastq1} -C fastq
            fastqs=~{find fastq -type f}
            readFilesCommand=""
            if [[ "$fastqs" = *.gz ]] ; then
                readFilesCommand="--readFilesCommand \"gunzip -c\""
            fi
        fi

    STAR \
        --runMode alignReads \
        --genomeDir $genomeDir \
        --genomeFastaFiles ~{insertions_fasta_file} \
        --runThreadN $cpu \
        --readFilesIn $fastqs \
        $readFilesCommand \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ~{sample_id}. \
        --outSAMstrandField intronMotif \
        --outSAMunmapped Within \
        ~{"--twopassMode " + two_pass_mode} \
        --alignSJDBoverhangMin 10 \
        --genomeSuffixLengthMax 10000 \
        --limitBAMsortRAM 47271261705 \
        --alignInsertionFlush Right \
        --outSAMfilter KeepOnlyAddedReferences \
        --alignMatesGapMax ${max_mate_dist} \
        --alignIntronMax  ${max_mate_dist} \
        --peOverlapNbasesMin 12 \
        --peOverlapMMp 0.1 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --alignSplicedMateMapLminOverLmate 0 \
        --alignSplicedMateMapLmin 30 \

        samtools index "${sample_id}.Aligned.sortedByCoord.out.bam"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
