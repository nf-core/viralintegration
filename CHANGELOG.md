# nf-core/viralintegration: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.1.1 - [2023-07-19]

### `Added`

- Pipeline summary to README
- Update to nf-core/tools v2.9

## v0.1.0 - [2023-03-29]

Initial release of nf-core/viralintegration, created with the [nf-core](https://nf-co.re/) template. (@alyssa-ab) (@Emiller88)

This pipeline is a re-implementation of [CTAT-VirusIntegrationFinder v1.5.0](https://github.com/broadinstitute/CTAT-VirusIntegrationFinder).

### `Pipeline Summary`

1. Input Check
   - Input path to sample FASTAs in samplesheet.csv
   - Check that sample meets requirements (samplesheet_check)
2. Read QC (FastQC)
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
9. Present quality checking and visualization for raw reads, adaptor trimming, and STAR alignments (MultiQC)

### `Added`

- Add CTAT-VIF virus_db.fasta [#1](https://github.com/nf-core/viralintegration/pull/1)
- Add small human test data set (chromosomes 6, 11, and 18 FASTA and GTF) [#31](https://github.com/nf-core/viralintegration/issues/31)
- Write nf-test for full workflow [#39](https://github.com/nf-core/viralintegration/issues/39)
- Add local module labels for resource management [#35](https://github.com/nf-core/viralintegration/issues/35)
- Write documentation [#50](https://github.com/nf-core/viralintegration/issues/50)
