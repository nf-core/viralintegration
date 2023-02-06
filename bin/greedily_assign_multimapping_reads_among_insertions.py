#!/usr/bin/env python

#########################################
# Author: [brianjohnhaas](https://github.com/brianjohnhaas)
# File: greedily_assign_multimapping_reads_among_insertions.py
# Source: https://github.com/broadinstitute/CTAT-VirusIntegrationFinder/blob/master/util/greedily_assign_multimapping_reads_among_insertions.py
# Source+commit: https://github.com/broadinstitute/CTAT-VirusIntegrationFinder/blob/551eff101679a0271d8bdfbdfb0d8724c2975c13/util/greedily_assign_multimapping_reads_among_insertions.py
# Download Date: 2022-12-28, commit: 551eff1
# This source code is licensed under the BSD 3-Clause license
#########################################

import sys, os, re
import csv
from collections import defaultdict
import argparse


def main():
    parser = argparse.ArgumentParser(
        description="regroups breakpoints based on virus single breakends",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--init_full_tsv",
        type=str,
        default="",
        required=True,
        help="intial (preliminary phase 1) insertion predictions",
    )
    parser.add_argument(
        "--MIN_ALT_BREAK_FRAC_READS",
        type=float,
        default=0.1,
        help="requires at least this fraction the alt breakpoint adj read count of the top scoring breakpoint for inclusion as alt.",
    )
    parser.add_argument(
        "--MAX_FRAC_MULTIMAPPING_NONPRIMARY",
        type=float,
        default=0.25,
        help="max fraction of multimapping reads allowed for non-primary insertion record",
    )
    parser.add_argument(
        "--include_readnames", action="store_true", default=False, help="include readnames for alignment evidence"
    )

    args = parser.parse_args()

    init_full_tsv = args.init_full_tsv
    MIN_ALT_BREAK_FRAC_READS = args.MIN_ALT_BREAK_FRAC_READS
    MAX_FRAC_MULTIMAPPING_NONPRIMARY = args.MAX_FRAC_MULTIMAPPING_NONPRIMARY
    INCLUDE_READNAMES = args.include_readnames

    # begin work
    csv.field_size_limit(sys.maxsize)
    tab_reader = csv.DictReader(open(init_full_tsv, "rt"), delimiter="\t")

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
    """

    rows = list()
    rows_grouped_by_virus_brkpt = defaultdict(list)
    for row in tab_reader:
        if not (
            (re.match("chr[\dMXY]+$", row["chrA"]) is not None) ^ (re.match("chr[\dMXY]+$", row["chrB"]) is not None)
        ):
            # only considering the main chromosomes
            continue

        row["total"] = int(row["total"])
        rows.append(row)
        virus, brkpt = (
            (row["chrA"], row["coordA"]) if re.match("chr[\dMXY]+$", row["chrB"]) else (row["chrB"], row["coordB"])
        )
        virus_brkpt_token = ":".join([virus, brkpt])
        rows_grouped_by_virus_brkpt[virus_brkpt_token].append(row)

    rows = sorted(rows, key=lambda x: (x["total"], x["entry"]), reverse=True)

    fieldnames = list(tab_reader.fieldnames) + [
        "adj_total",
        "excluded_reads",
        "frac_reads_removed",
        "virus_brkend_grp",
        "is_primary",
    ]
    if not INCLUDE_READNAMES:
        fieldnames.remove("readnames")
        fieldnames.remove("excluded_reads")

    writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()

    def sort_ranking(row):
        total = row["total"]
        is_HPV = 1 if re.search("HPV", row["entry"]) else 0
        # in case of ties, prefer the HPV accession entry over an alternative representation
        return (total, is_HPV)

    processed_brkpts = set()
    ev_reads_seen = set()

    rows = sorted(rows, key=sort_ranking, reverse=True)

    for row in rows:
        virus, brkpt = (
            (row["chrA"], row["coordA"]) if re.match("chr[\dMXY]+$", row["chrB"]) else (row["chrB"], row["coordB"])
        )
        virus_brkpt_token = ":".join([virus, brkpt])
        if virus_brkpt_token in processed_brkpts:
            continue

        shared_virus_brkpt_rows = rows_grouped_by_virus_brkpt[virus_brkpt_token]

        # define adjusted totals based on earlier defined multimapping reads
        for loc_row in shared_virus_brkpt_rows:
            compute_adjusted_total(loc_row, ev_reads_seen)

        shared_virus_brkpt_rows = sorted(shared_virus_brkpt_rows, key=lambda x: x["adj_total"], reverse=True)

        top_scoring_brkpt_row = shared_virus_brkpt_rows[0]
        top_score = top_scoring_brkpt_row["adj_total"]
        top_scoring_brkpt_row["virus_brkend_grp"] = virus_brkpt_token
        top_scoring_brkpt_row["is_primary"] = "True"

        add_ev_read_exclusion(top_scoring_brkpt_row, ev_reads_seen)

        for remaining_row in shared_virus_brkpt_rows[1:]:
            add_ev_read_exclusion(remaining_row, ev_reads_seen)
            remaining_row["virus_brkend_grp"] = virus_brkpt_token
            remaining_row["is_primary"] = "False"

        # report entries
        for loc_row in shared_virus_brkpt_rows:
            if not INCLUDE_READNAMES:
                del loc_row["readnames"]
                del loc_row["excluded_reads"]

            if (
                loc_row["adj_total"] > 0
                and loc_row["adj_total"] / top_score >= MIN_ALT_BREAK_FRAC_READS
                and float(loc_row["frac_reads_removed"]) < MAX_FRAC_MULTIMAPPING_NONPRIMARY
            ):
                writer.writerow(loc_row)

        processed_brkpts.add(virus_brkpt_token)


def compute_adjusted_total(row, ev_reads_seen):
    readnames = row["readnames"].split(",")
    adj_readnames = list()
    excluded_readnames = list()
    for readname in readnames:
        if readname in ev_reads_seen:
            excluded_readnames.append(readname)
        else:
            adj_readnames.append(readname)

    row["adj_total"] = len(adj_readnames)
    row["readnames"] = ",".join(adj_readnames)
    row["excluded_reads"] = ",".join(excluded_readnames)
    row["frac_reads_removed"] = "{:.3f}".format(len(excluded_readnames) / len(readnames))


def add_ev_read_exclusion(row, ev_reads_seen):
    readnames = row["readnames"].split(",")
    for readname in readnames:
        ev_reads_seen.add(readname)


if __name__ == "__main__":
    main()
