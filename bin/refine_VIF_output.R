#!/usr/bin/env Rscript

#########################################
# Author: [brianjohnhaas](https://github.com/brianjohnhaas)
# File: refine_VIF_output.R
# Source: https://github.com/broadinstitute/CTAT-VirusIntegrationFinder/blob/master/util/refine_VIF_output.Rscript
# Source+commit: https://github.com/broadinstitute/CTAT-VirusIntegrationFinder/blob/33dcbb06b8394e973aa01c536363557689cea23c/util/refine_VIF_output.Rscript
# Download Date: 2022-12-28, commit: 33dcbb0
# This source code is licensed under the BSD 3-Clause license
#########################################


suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("tidyverse"))

parser <- ArgumentParser()
parser$add_argument("--prelim_counts", help = "prelim counts file w stats", required = TRUE, nargs = 1)
parser$add_argument("--vif_counts", help = "vif counts file", required = TRUE)
parser$add_argument("--output", help = "output filename", required = TRUE)

args <- parser$parse_args()


prelim_counts_filename <- args$prelim_counts
vif_counts_filename <- args$vif_counts
output_filename <- args$output


prelim_counts <- read.table(prelim_counts_filename, header = T, sep = "\t", stringsAsFactors = F, comment.char = "", quote = "")
prelim_counts <- prelim_counts %>% rename(
    contig = entry,
    prelim.total = total,
    prelim.adj_total = adj_total,
    prelim.frac_reads_removed = frac_reads_removed,
    prelim.primary_brkpt_type = primary_brkpt_type,
    prelim.num_primary_reads = num_primary_reads,
    prelim.num_supp_reads = num_supp_reads,
    prelim.mean_hits = hits,
    prelim.mean_min_per_id = min_per_id,
    prelim.mean_max_end_clipping = max_end_clipping,
    prelim.mean_min_anchor_len = min_anchor_len
)

prelim_counts <- prelim_counts %>% select(
    contig, chrA, coordA, orientA, chrB, coordB, orientB,
    prelim.primary_brkpt_type, prelim.num_primary_reads, prelim.num_supp_reads,
    prelim.mean_hits, prelim.mean_min_per_id, prelim.mean_max_end_clipping, prelim.mean_min_anchor_len,
    prelim.total, prelim.adj_total, prelim.frac_reads_removed,
    entropyA, entropyB, flankA, flankB, splice_type,
    virus_brkend_grp, is_primary
)

vif_counts <- read.table(vif_counts_filename, header = T, sep = "\t", stringsAsFactors = F, comment.char = "", quote = "")

data <- right_join(prelim_counts, vif_counts, by = "contig")

# data[is.na(data)] = 0

# join with non-primary insetions of the same virus brkpt group
data <- data %>% mutate(search_grp = paste(virus_brkend_grp, flankA, flankB))
virus_brkend_grps_want <- data %>% pull(search_grp)
prelim_counts <- prelim_counts %>% mutate(search_grp = paste(virus_brkend_grp, flankA, flankB))

# join with imputed entries
data <- bind_rows(
    data,
    prelim_counts %>% filter(is_primary == "False") %>% filter(search_grp %in% virus_brkend_grps_want)
)

# for the imputed entries, assign total read count to that of the virus brken grp
data <- data %>%
    group_by(search_grp) %>%
    mutate(total2 = max(total, na.rm = T)) %>%
    ungroup()
data <- data %>%
    mutate(total = ifelse(is.na(total), total2, total)) %>%
    select(-total2)

data <- data %>% select(-search_grp)

data <- data %>% arrange(desc(total), desc(prelim.adj_total))

# data = data %>% arrange(search_grp)

write.table(data, file = output_filename, row.names = F, sep = "\t", quote = F)


quit(save = "no", status = 0, runLast = FALSE)
