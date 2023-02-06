#!/usr/bin/env Rscript

#########################################
# Author: [brianjohnhaas](https://github.com/brianjohnhaas)
# File: make_VIF_genome_abundance_plot.R
# Source: https://github.com/broadinstitute/CTAT-VirusIntegrationFinder/blob/master/util/make_VIF_genome_abundance_plot.Rscript
# Source+commit: https://github.com/broadinstitute/CTAT-VirusIntegrationFinder/blob/4a5502f1c1421e4f2d3eb2295da5e970fce46e13/util/make_VIF_genome_abundance_plot.Rscript
# Download Date: 2022-12-28, commit: 4a5502f
# This source code is licensed under the BSD 3-Clause license
#########################################

options(bitmapType = "cairo") # needed to avoid X11 issue

condense_data <- function(data) {
    data_condensed <- bind_rows(
        data %>% select(chr = chrA, coord = coordA, orient = orientA, counts = total),
        data %>% select(chr = chrB, coord = coordB, orient = orientB, counts = total)
    )

    data_condensed$chr <- str_replace(data_condensed$chr, "chr", "")

    return(data_condensed)
}


chroms <- c(c(1:22), "X", "Y", "M")

chr_lengths <- data.frame(
    chr = factor(chroms, levels = chroms),
    chr_begin = 0,
    chr_length = c(
        248956422,
        242193529,
        198295559,
        190214555,
        181538259,
        170805979,
        159345973,
        145138636,
        138394717,
        133797422,
        135086622,
        133275309,
        114364328,
        107043718,
        101991189,
        90338345,
        83257441,
        80373285,
        58617616,
        64444167,
        46709983,
        50818468,
        156040895,
        57227415,
        16569
    )
)

get_chr_panel_base_plot <- function(dataToPlot, chr_lens = chr_lengths) {
    dataToPlot <- dataToPlot %>% filter(chr %in% chroms)

    dataToPlot$chr <- factor(dataToPlot$chr, levels = chroms)

    p <- ggplot(data = dataToPlot) +
        facet_grid(~chr, scales = "free_x", space = "fixed") +
        theme_bw() +
        theme(
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank()
        ) +
        geom_vline(data = chr_lens, aes(xintercept = chr_begin), color = NA) +
        geom_vline(data = chr_lens, aes(xintercept = chr_length), color = NA)

    return(p)
}



make_chr_plot <- function(dataToPlot, title) {
    p <- get_chr_panel_base_plot(dataToPlot)

    p <- p + geom_point(aes(x = coord, y = counts)) + ggtitle(title)

    return(p)
}



main <- function() {
    suppressPackageStartupMessages(library("argparse"))
    suppressPackageStartupMessages(library("tidyverse"))

    parser <- ArgumentParser()
    parser$add_argument("--vif_report", help = "final counts file", required = TRUE, nargs = 1)
    parser$add_argument("--output_png", help = "output png filename", required = TRUE)
    parser$add_argument("--title", help = "output png title", required = TRUE)

    args <- parser$parse_args()
    vif_report <- args$vif_report
    title <- args$title
    output_png <- args$output_png

    data <- read.table(vif_report, header = T, sep = "\t", stringsAsFactors = F, comment.char = "", quote = "")
    data <- condense_data(data)
    p <- make_chr_plot(data, title)
    png(output_png, width = 800, height = 400)
    plot(p)
    dev.off()

    quit(save = "no", status = 0, runLast = FALSE)
}


if (length(sys.calls()) == 0) {
    main()
}
