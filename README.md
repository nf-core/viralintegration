# ![nf-core/viralintegration](docs/images/nf-core-viralintegration_logo_light.png#gh-light-mode-only) ![nf-core/viralintegration](docs/images/nf-core-viralintegration_logo_dark.png#gh-dark-mode-only)

[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/viralintegration/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.7783480-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.7783480)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/nf-core/viralintegration)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23viralintegration-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/viralintegration)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/viralintegration** is a bioinformatics best-practice analysis pipeline for the identification of viral integration events in genomes using a chimeric read approach. It was initially based on the [CTAT-VirusIntegrationFinder](https://github.com/broadinstitute/CTAT-VirusIntegrationFinder).

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

<!-- TODO nf-core: Add full-sized test dataset and amend the paragraph below if applicable -->

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources.The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/viralintegration/results).

## Pipeline summary

1. Input Check
   - Input path to sample FASTAs in samplesheet.csv
   - Check that sample meets requirements (samplesheet_check)
2. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
3. Align reads to human genome
   - Generate index and perform alignment (STAR)
4. Quality trimming for unaligned reads
   - Quality and adaptor trimming (Trimmomatic)
   - Remove polyAs from reads (PolyAStripper)
5. Identify chimeric reads
   - Combine human and virus FASTAs (cat_fasta)
   - Generate index and perform alignment to combined human + viral reference (STAR)
   - Sort and index alignments (SAMtools)
   - Determine potential insertion site candidates and optimize file (insertion_site_candidates, abridged_TSV)
6. Virus Report outputs:
   - Viral read counts in a tsv table and png plot
   - Preliminary genome wide abundance plot
   - Bam and bai for reads detected in potential viral insertion site
   - Web based interactive genome viewer for virus infection evidence (VirusDetect.igvjs.html)
7. Verify chimeric reads
   - Create chimeric FASTA and GTF extracts (extract_chimeric_genomic_targets)
   - Generate index and perform alignment to verify chimeric reads (STAR)
   - Sort and index validated alignments (SAMtools)
   - Remove duplicate alignments (remove_duplicates)
   - Generate evidence counts for chimeric reads (chimeric_contig_evidence_analyzer)
8. Summary Report outputs:
   - Refined genome wide abundance plog png
   - Insertion site candidates in tab-delimited format with gene annotations (vif.refined.wRefGeneAnnots.tsv)
   - Web based interactive genome viewer for virus insertion sites (vif.html)
9. Present quality checking and visualization for raw reads, adaptor trimming, and STAR alignments ([`MultiQC`](http://multiqc.info/))

## Usage

> **Note**
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how
> to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline)
> with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

Now, you can run the pipeline using:

```bash
nextflow run nf-core/viralintegration \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> **Warning:**
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those
> provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/viralintegration/usage) and the [parameter documentation](https://nf-co.re/viralintegration/parameters).

```bash
nextflow run nf-core/viralintegration --input samplesheet.csv --outdir <OUTDIR> --genome GRCh37 -profile <docker/singularity/podman/shifter/charliecloud/conda/institute>
```

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/viralintegration/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/viralintegration/output).

## Credits

nf-core/viralintegration was originally written by Alyssa Briggs ([@alyssa-ab](https://github.com/alyssa-ab)) and Edmund Miller ([@Emiller88](https://github.com/emiller88)) from [The Functional Genomics Laboratory](https://taehoonkim.org/) at [The Univeristy of Texas at Dallas](https://www.utdallas.edu/).

We thank the following people for their extensive assistance in the development of this pipeline:

- [Tae Hoon Kim](https://github.com/taehoonkim-phd)
- [Robert Allaway](https://github.com/allaway)

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#viralintegration` channel](https://nfcore.slack.com/channels/viralintegration) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

If you use nf-core/viralintegration for your analysis, please cite it using the following doi: [10.5281/zenodo.7783480](https://doi.org/10.5281/zenodo.7783480)

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
