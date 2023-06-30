#!/usr/bin/env python

#########################################
# Author: [brianjohnhaas](https://github.com/brianjohnhaas)
# File: extract_prelim_chimeric_genome_read_alignments.py
# Source: https://github.com/broadinstitute/CTAT-VirusIntegrationFinder/blob/master/util/extract_prelim_chimeric_genome_read_alignments.py
# Source+commit: https://github.com/broadinstitute/CTAT-VirusIntegrationFinder/blob/1553104ae8fb5bd1631d8122ff1b5cc9e6093640/util/extract_prelim_chimeric_genome_read_alignments.py
# Download Date: 2022-12-28, commit: 1553104
# This source code is licensed under the BSD 3-Clause license
#########################################

import sys, os, re
import pysam
from collections import defaultdict
import logging
import argparse
import pandas as pd
import subprocess


logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description="add alignment stats", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("--star_bam", required=True, type=str, help="STAR read alignments to genome+virus db")
    parser.add_argument("--vif_full_tsv", required=True, type=str, help="VIF full tsv containing evidence read names")
    parser.add_argument(
        "--output_bam", required=True, type=str, help="output bam containing just the chimeric read alignments"
    )

    args = parser.parse_args()

    bam_filename = args.star_bam
    vif_full_tsv = args.vif_full_tsv
    output_bam_filename = args.output_bam

    samfile = pysam.AlignmentFile(bam_filename, "rb")

    bam_outfile = pysam.AlignmentFile(output_bam_filename, "wb", template=samfile)

    ev_reads = set()

    not_seen = ev_reads.copy()

    logger.info("-capturing reads of interest from {}".format(vif_full_tsv))
    vif_df = pd.read_csv(vif_full_tsv, sep="\t")
    vif_df = vif_df.astype({"readnames": "str"})  # in case they look like integers.
    for _, row in vif_df.iterrows():
        readnames = row["readnames"].split(",")
        for readname in readnames:
            ev_reads.add(readname)

    logger.info("-capturing read alignments from {}, writing to: {}".format(bam_filename, output_bam_filename))

    read_counter = 0

    for read in samfile.fetch():
        read_name = read.query_name

        if read_name not in ev_reads:
            continue

        read_counter += 1

        if read_counter % 100 == 0:
            logger.info("-captured {} chimeric alignments".format(read_counter))

        bam_outfile.write(read)

        if read_name in not_seen:
            not_seen.remove(read_name)

    if not_seen:
        raise RuntimeError(
            "Error, missing alignments for {} chimeric reads: {}".format(len(not_seen), ",".join(list(not_seen)))
        )

    bam_outfile.close()
    samfile.close()

    subprocess.check_call("samtools index {}".format(output_bam_filename), shell=True)

    logger.info("-done")

    sys.exit(0)


if __name__ == "__main__":
    main()
