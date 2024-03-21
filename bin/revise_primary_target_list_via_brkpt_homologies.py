#!/usr/bin/env python

#########################################
# Author: [brianjohnhaas](https://github.com/brianjohnhaas)
# File: revise_primary_target_list_via_brkpt_homologies.py
# Source: https://github.com/broadinstitute/CTAT-VirusIntegrationFinder/blob/master/util/revise_primary_target_list_via_brkpt_homologies.py
# Source+commit: https://github.com/broadinstitute/CTAT-VirusIntegrationFinder/blob/e20a1eb2acefeb23e36719f4c03b0e4e93e6a754/util/revise_primary_target_list_via_brkpt_homologies.py
# Download Date: 2022-12-28, commit: e20a1eb
# This source code is licensed under the BSD 3-Clause license
#########################################

import sys, os, re
import csv
from collections import defaultdict
import argparse


def main():
    parser = argparse.ArgumentParser(
        description="assign Maybe status to top entries of breakpoint homology groups",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--vif_tsv", type=str, default="", required=True, help="intial (preliminary phase 1) insertion predictions"
    )

    args = parser.parse_args()

    vif_tsv = args.vif_tsv

    # begin work
    csv.field_size_limit(sys.maxsize)
    tab_reader = csv.DictReader(open(vif_tsv, "rt"), delimiter="\t")

    """
    0       entry
    1       chrA
    2       coordA
    3       orientA
    4       chrB
    5       coordB
    6       orientB
    7       primary_brkpt_type
    8       num_primary_reads
    9       num_supp_reads
    10      total
    11      readnames
    12      adj_total
    13      excluded_reads
    14      frac_reads_removed
    15      virus_brkend_grp
    16      is_primary
    17      hits
    18      min_per_id
    19      max_end_clipping
    20      min_anchor_len
    21      flankA
    22      flankB
    23      entropyA
    24      entropyB
    25      splice_type
    """

    rows = list()
    rows_grouped_by_virus_brkpt = defaultdict(list)
    for row in tab_reader:
        row["total"] = int(row["total"])
        rows.append(row)
        virus, brkpt = (
            (row["chrA"], row["coordA"]) if re.match("chr[\dMXY]+$", row["chrB"]) else (row["chrB"], row["coordB"])
        )
        virus_brkpt_token = ":".join([virus, brkpt, row["flankA"], row["flankB"]])
        rows_grouped_by_virus_brkpt[virus_brkpt_token].append(row)

    rows = sorted(rows, key=lambda x: (x["total"], x["entry"]), reverse=True)

    fieldnames = list(tab_reader.fieldnames)
    writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()

    processed_brkpts = set()
    ev_reads_seen = set()

    for row in rows:
        virus, brkpt = (
            (row["chrA"], row["coordA"]) if re.match("chr[\dMXY]+$", row["chrB"]) else (row["chrB"], row["coordB"])
        )
        virus_brkpt_token = ":".join([virus, brkpt, row["flankA"], row["flankB"]])

        if virus_brkpt_token in processed_brkpts:
            continue

        shared_virus_brkpt_rows = rows_grouped_by_virus_brkpt[virus_brkpt_token]

        shared_virus_brkpt_rows = sorted(shared_virus_brkpt_rows, key=lambda x: x["adj_total"], reverse=True)

        if shared_virus_brkpt_rows[0]["is_primary"] == "False":
            # pursue it as an alt candidate
            shared_virus_brkpt_rows[0]["is_primary"] = "Maybe"

        # report entries
        for loc_row in shared_virus_brkpt_rows:
            writer.writerow(loc_row)

        processed_brkpts.add(virus_brkpt_token)


if __name__ == "__main__":
    main()
