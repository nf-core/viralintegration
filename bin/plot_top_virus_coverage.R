#!/usr/bin/env Rscript

#########################################
# Author: [brianjohnhaas](https://github.com/brianjohnhaas)
# File: plot_top_virus_coverage.R
# Source: https://github.com/broadinstitute/CTAT-VirusIntegrationFinder/blob/master/util/plot_top_virus_coverage.Rscript
# Source+commit: https://github.com/broadinstitute/CTAT-VirusIntegrationFinder/blob/caddb8a67c6de0b9c3ca4b61d1a9833d94f619c7/util/plot_top_virus_coverage.Rscript
# Download Date: 2022-12-28, commit: caddb8a
# This source code is licensed under the BSD 3-Clause license
#########################################

options(bitmapType = "cairo") # needed to avoid X11 issue

suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("tidyverse"))
parser <- ArgumentParser()
parser$add_argument("--vif_report", help = "final counts file", required = TRUE, nargs = 1)
parser$add_argument("--bam", help = "alignments in sorted bam file", required = TRUE, nargs = 1)
parser$add_argument("--virus_fai", help = "virus fasta fai file", required = TRUE, nargs = 1)
parser$add_argument("--output_prefix", help = "output prefix for png and tsv filenames", required = TRUE, nargs = 1)
parser$add_argument("--depth", action = "store_true", dest = "depth", help = "Create depth plots")

args <- parser$parse_args()

vif_report <- args$vif_report
virus_seqlens_file <- args$virus_seqlens
bam_file <- args$bam
virus_fai_file <- args$virus_fai
output_prefix <- args$output_prefix
depth_plots <- args$depth

library(tidyverse)


run_command <- function(cmd) {
    message(cmd)
    ret <- system(cmd)
    if (ret != 0) {
        stop(paste("Error, cmd: ", cmd, " died wit ret ", ret))
    }
}




bam_index <- paste0(bam_file, ".bai")
if (!file_test("-f", bam_index)) {
    ## must create it.
    cmd <- paste("samtools index", bam_file)
    run_command(cmd)
}


viruses <- read.table(virus_fai_file, header = F, stringsAsFactors = F, sep = "\t")[, 1]


get_virus_mappings <- function(data) {
    data_condensed <- bind_rows(
        data %>% select(chr = chrA, coord = coordA, orient = orientA, counts = total),
        data %>% select(chr = chrB, coord = coordB, orient = orientB, counts = total)
    )

    # data_condensed$chr = str_replace(data_condensed$chr, "chr", "")

    virus_mappings <- data_condensed %>% filter(chr %in% viruses)

    return(virus_mappings)
}


data <- read.table(vif_report, header = T, stringsAsFactors = FALSE, sep = "\t", comment.char = "", quote = "")
counts_per_virus <- NULL
if (nrow(data) > 0) {
    virus_mappings <- get_virus_mappings(data)
    counts_per_virus <- virus_mappings %>%
        rename(virus = chr) %>%
        group_by(virus) %>%
        summarize_at(vars(counts), list(chim_reads = sum))
}

idxstats_filename <- paste0(bam_file, ".idxstats")
if (!file_test("-f", idxstats_filename)) {
    cmd <- paste("samtools idxstats", bam_file, ">", idxstats_filename)
    run_command(cmd)
}

idxstats <- read.table(idxstats_filename, header = F, sep = "\t", stringsAsFactors = F, comment.char = "", quote = "")
colnames(idxstats) <- c("virus", "seqlen", "mapped", "unmapped")
idxstats <- idxstats %>% filter(virus %in% viruses)


if (!is.null(counts_per_virus)) {
    virus_count_info <- full_join(idxstats, counts_per_virus, by = "virus")
} else {
    virus_count_info <- idxstats %>% mutate(chim_reads = 0)
}

virus_count_info <- virus_count_info %>%
    filter(mapped > 0 | chim_reads > 0) %>%
    select(virus, seqlen, mapped, chim_reads) %>%
    arrange(desc(mapped))
virus_count_info[is.na(virus_count_info)] <- 0

png(paste0(output_prefix, ".virus_read_counts.png"), width = 800, height = 400)
p <- virus_count_info %>%
    gather(key = map_type, value = read_count, c(mapped, chim_reads)) %>%
    ggplot(aes(x = virus, y = read_count)) +
    geom_col(aes(fill = map_type), position = "dodge") +
    ggtitle("Virus-mapped and chimeric read counts")
plot(p)
dev.off()


png(paste0(output_prefix, ".virus_read_counts_log.png"), width = 800, height = 400)
p <- virus_count_info %>%
    gather(key = map_type, value = read_count, c(mapped, chim_reads)) %>%
    ggplot(aes(x = virus, y = read_count)) +
    geom_col(aes(fill = map_type), position = "dodge") +
    scale_y_log10() +
    ggtitle("Virus-mapped and chimeric read counts (log10scale)")
plot(p)
dev.off()



virus_count_info <- virus_count_info %>% filter(mapped > 0)
if (nrow(virus_count_info) == 0) {
    message("no virus coverage info to report")
    quit(save = "no", status = 0, runLast = FALSE)
}




read_depth_df <- NULL

for (virus in virus_count_info$virus) {
    message("-examining coverage depth for virus:", virus)
    depth_file <- "depth.tsv"
    depth_script_py <- "sam_depth_ignore_gaps.py"
    cmd <- paste("bash -c \'set -eou pipefail && samtools view -h ", bam_file, paste0('"', virus, '"'), " | ", depth_script_py, " - >", depth_file, "'")
    run_command(cmd)
    n_bases_covered <- 0

    if (file.size(depth_file) > 0) {
        depth_data <- read.table(depth_file, header = F, sep = "\t", comment.char = "", quote = "")[, c(2, 3)]
        n_bases_covered <- nrow(depth_data)
        colnames(depth_data) <- c("pos", "depth")

        seqlen <- idxstats$seqlen[which(idxstats$virus == virus)]

        allpos <- data.frame(pos = 1:seqlen)

        depth_data <- left_join(allpos, depth_data, by = "pos")
        depth_data$depth[is.na(depth_data$depth)] <- 0

        if (depth_plots) {
            p <- depth_data %>% ggplot(aes(x = pos, y = depth)) +
                geom_line() +
                ggtitle(paste0("Viral Genome Depth of Coverage for [", virus, "]"))
            png(paste0(output_prefix, ".virus_coverage_", virus, ".png"), width = 800, height = 400)
            plot(p)
            dev.off()
        }
    }
    read_depth_df <- bind_rows(read_depth_df, data.frame(virus = virus, n_bases_covered = n_bases_covered))
}


virus_count_info <- left_join(virus_count_info, read_depth_df, by = "virus")
virus_count_info <- virus_count_info %>%
    mutate(frac_covered = n_bases_covered / seqlen) %>%
    filter(frac_covered > 1e-3)
virus_count_info <- virus_count_info %>%
    arrange(desc(frac_covered)) %>%
    mutate(frac_covered = sprintf("%.3f", frac_covered))


virus_mappings_summary_file <- paste0(output_prefix, ".virus_read_counts_summary.tsv")
write.table(virus_count_info, file = virus_mappings_summary_file, quote = F, sep = "\t", row.names = F)

quit(save = "no", status = 0, runLast = FALSE)
