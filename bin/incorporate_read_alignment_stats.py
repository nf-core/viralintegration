#!/usr/bin/env python

#########################################
# Author: [brianjohnhaas](https://github.com/brianjohnhaas)
# File: incorporate_read_alignment_stats.py
# Source: https://github.com/broadinstitute/CTAT-VirusIntegrationFinder/blob/master/util/incorporate_read_alignment_stats.py
# Source+commit: https://github.com/broadinstitute/CTAT-VirusIntegrationFinder/blob/9af47f4567cca1263404aa0c535125f34a7577cc/util/incorporate_read_alignment_stats.py
# Download Date: 2022-12-28, commit: 9af47f4
# This source code is licensed under the BSD 3-Clause license
#########################################

import sys, os, re
import pysam
from collections import defaultdict
import logging
import argparse
import pandas as pd
import statistics as st

logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description="add alignment stats", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("--supp_reads_bam", required=True, type=str, help="supplemental read alignments")
    parser.add_argument("--vif_full_tsv", required=True, type=str, help="VIF full tsv containing evidence read names")
    parser.add_argument("--output", required=True, type=str, help="output tsv file containing read alignment stats")
    parser.add_argument(
        "--detailed", action="store_true", help="show attribute values for each read instead of mean statistic"
    )

    args = parser.parse_args()

    bam_filename = args.supp_reads_bam
    vif_full_tsv = args.vif_full_tsv
    outputfilename = args.output

    samfile = pysam.AlignmentFile(bam_filename, "rb")

    ev_reads = set()

    logger.info("-capturing reads of interest from {}".format(vif_full_tsv))
    vif_df = pd.read_csv(vif_full_tsv, sep="\t")
    vif_df = vif_df.astype({"readnames": "str"})  # in case they look like integers
    for _, row in vif_df.iterrows():
        readnames = row["readnames"].split(",")
        for readname in readnames:
            ev_reads.add(readname)

    logger.info("-capturing read alignment stats from bam file: {}".format(bam_filename))

    read_to_hit_count = defaultdict(int)
    read_to_min_per_id = dict()
    read_to_max_end_clipping = defaultdict(int)
    read_to_min_anchor_len = dict()

    read_counter = 0

    for read in samfile.fetch():
        read_name = read.query_name

        if read_name not in ev_reads:
            continue

        read_counter += 1

        if read_counter % 10000 == 0:
            logger.info("-processed {} alignments".format(read_counter))

        aligned_bases = len(read.get_aligned_pairs(matches_only=True))
        if read_name in read_to_min_anchor_len:
            read_to_min_anchor_len[read_name] = min(aligned_bases, read_to_min_anchor_len[read_name])
        else:
            read_to_min_anchor_len[read_name] = aligned_bases

        NH = read.get_tag("NH")
        read_to_hit_count[read_name] = max(read_to_hit_count[read_name], NH)

        mismatch_count = read.get_tag("NM")
        per_id = 100.0 - (float(mismatch_count) / aligned_bases * 100.0)

        if read_name in read_to_min_per_id:
            read_to_min_per_id[read_name] = min(read_to_min_per_id[read_name], per_id)
        else:
            read_to_min_per_id[read_name] = per_id

        cigar = read.cigarstring

        max_clip = 0
        m = re.search("^(\d+)S", cigar)
        if m:
            max_clip = int(m.group(1))

        m = re.search("(\d+)S$", cigar)
        if m:
            max_clip = max(max_clip, int(m.group(1)))

        read_to_max_end_clipping[read_name] = max(read_to_max_end_clipping[read_name], max_clip)

    logger.info("-generating alignment stats report")

    vif_df["hits"] = ""
    vif_df["min_per_id"] = ""
    vif_df["max_end_clipping"] = ""
    vif_df["min_anchor_len"] = ""

    for i, row in vif_df.iterrows():
        readnames = row["readnames"].split(",")
        hits = list()
        min_per_ids = list()
        max_end_clipping = list()
        min_anchor_lengths = list()

        for readname in readnames:
            if readname not in read_to_hit_count:
                raise RuntimeError("Error, missing hit count for read: {}".format(readname))

            hits.append(read_to_hit_count[readname])
            min_per_ids.append(read_to_min_per_id[readname])
            max_end_clipping.append(read_to_max_end_clipping[readname])
            min_anchor_lengths.append(read_to_min_anchor_len[readname])

        if args.detailed:
            vif_df.loc[i, "hits"] = ",".join([str(x) for x in hits])
            vif_df.loc[i, "min_per_id"] = ",".join(["{:.1f}".format(x) for x in min_per_ids])
            vif_df.loc[i, "max_end_clipping"] = ",".join([str(x) for x in max_end_clipping])
            vif_df.loc[i, "min_anchor_len"] = ",".join([str(x) for x in min_anchor_lengths])

        else:
            vif_df.loc[i, "hits"] = "{:.3f}".format(st.mean(hits))
            vif_df.loc[i, "min_per_id"] = "{:.1f}".format(st.mean(min_per_ids))
            vif_df.loc[i, "max_end_clipping"] = "{:.3f}".format(st.mean(max_end_clipping))
            vif_df.loc[i, "min_anchor_len"] = "{:.3f}".format(st.mean(min_anchor_lengths))

    # vif_df.drop('readnames', axis=1, inplace=True)

    logger.info("-writing outputfile: {}".format(outputfilename))
    vif_df.to_csv(outputfilename, sep="\t", index=False)

    logger.info("-done")

    sys.exit(0)


if __name__ == "__main__":
    main()
