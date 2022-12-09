#!/usr/bin/env python

import sys, os, re
import csv
from collections import defaultdict
import argparse



def main():


    parser = argparse.ArgumentParser(description="assign Maybe status to top entries of breakpoint homology groups", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--vif_tsv", type=str, default="", required=True, help="final prediction list including primary and non-primary entries")

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

    rows_grouped_by_virus_brkpt = defaultdict(list)

    for row in tab_reader:
        virus_brkpt_token = ":".join([row['virus_brkend_grp'], row['flankA'], row['flankB'] ])
        rows_grouped_by_virus_brkpt[virus_brkpt_token].append(row)


    final_rows = list()

    for virus_brkpt_token in rows_grouped_by_virus_brkpt:
        rows = rows_grouped_by_virus_brkpt[virus_brkpt_token]

        primary_row = None
        nonprimary_info = list()

        for row in rows:
            if row['is_primary'] != "False":
                primary_row = row
            else:
                nonprimary_info.append(":".join([row['contig'], row['prelim.adj_total'] ]))

        assert primary_row is not None, "couldn't locate primary insertion entry for {}".format(virus_brkpt_token)

        primary_row['num_alt_locs'] = len(nonprimary_info)
        primary_row['alt_locs'] = ",".join(nonprimary_info)

        final_rows.append(primary_row)


    final_rows = sorted(final_rows, key=lambda x: x['total'], reverse=True)


    fieldnames = list(tab_reader.fieldnames) + ['num_alt_locs', 'alt_locs']
    fieldnames.remove('is_primary')
    writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()


    for row in final_rows:
        del row['is_primary']
        writer.writerow(row)


if __name__=='__main__':
    main()
