#!/usr/bin/env python3

#########################################
# Author: [brianjohnhaas](https://github.com/brianjohnhaas)
# File: pre_filter_non_human_virus_chimeric_alignments.py
# Source: https://github.com/broadinstitute/CTAT-VirusIntegrationFinder/blob/master/util/pre_filter_non_human_virus_chimeric_alignments.py
# Source+commit: https://github.com/broadinstitute/CTAT-VirusIntegrationFinder/blob/3938c9074bfe4a718df3f4067b328ab9e304b5a2/util/pre_filter_non_human_virus_chimeric_alignments.py
# Download Date: 2022-12-28, commit: 3938c90
# This source code is licensed under the BSD 3-Clause license
#########################################

import sys, os, re
import argparse
import subprocess
import logging
from collections import defaultdict
import pandas as pd
import csv
import gzip

if sys.version_info[0] < 3:
    raise Exception("This script requires Python 3")

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


"""  Chimeric.out.junction file formatting from the STAR manual:
        # from star doc:
        #The rst 9 columns give information about the chimeric junction:
        #The format of this le is as follows. Every line contains one chimerically aligned read, e.g.:
        #   chr22 23632601 + chr9 133729450 + 1 0 0 SINATRA-0006:3:3:6387:56650 23632554 47M29S 133729451 47S29M40p76M

        #The first 9 columns give information about the chimeric junction:
        #column 1: chromosome of the donor
        #column 2: rst base of the intron of the donor (1-based)
        #column 3: strand of the donor
        #column 4: chromosome of the acceptor
        #column 5: rst base of the intron of the acceptor (1-based)
        #column 6: strand of the acceptor
        #column 7: junction type: -1=encompassing junction (between the mates), 1=GT/AG, 2=CT/AC
        #column 8: repeat length to the left of the junction
        #column 9: repeat length to the right of the junction
        #Columns 10-14 describe the alignments of the two chimeric segments, it is SAM like. Alignments are given with respect to the (+)
 strand
        #column 10: read name
        #column 11: rst base of the rst segment (on the + strand)
        #column 12: CIGAR of the rst segment
        #column 13: rst base of the second segment
        #column 14: CIGAR of the second segment
"""


def main():
    parser = argparse.ArgumentParser(
        description="filters out chimeric alignments not involving human--virus or any involving chrM",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("--chimJ", type=str, required=True, help="STAR Chimeric.out.junction file")

    parser.add_argument(
        "--viral_db_fasta",
        type=str,
        required=True,
        help="database of the additional viral genome targets",
    )

    parser.add_argument(
        "--output",
        "-o",
        dest="output_filename",
        type=str,
        required=True,
        help="output filename",
    )
    parser.add_argument("--debug", action="store_true", help="debug mode")

    ###########################
    # Parse input arguments
    ###########################
    args_parsed = parser.parse_args()
    # Create constants
    chimJ_filename = args_parsed.chimJ
    viral_db_fasta_filename = args_parsed.viral_db_fasta
    output_filename = args_parsed.output_filename

    if args_parsed.debug:
        logger.setLevel(logging.DEBUG)

    ###########################
    ## get list of viral_db entries. (Viruses names found in fasta file)
    ###########################
    viral_db_entries = set()
    # cache the regular expression search, faster because more memory conservative
    match_regex = re.compile("^>(\S+)")
    with open(viral_db_fasta_filename, "rt") as fh:
        for line in fh:
            m = match_regex.search(line)
            if m:
                acc = m.group(1)
                viral_db_entries.add(acc)

    if os.path.splitext(chimJ_filename)[1] == ".gz":
        fh = gzip.open(chimJ_filename, "rt")
    else:
        fh = open(chimJ_filename, "rt")

    reader = csv.DictReader(fh, delimiter="\t")
    fieldnames = list(reader.fieldnames)

    ofh = open(output_filename, "wt")
    writer = csv.DictWriter(ofh, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()

    count_virus_chims = 0

    prev_read_name = ""
    reads_collected = list()
    prev_nonrelevant_chim = False

    for row in reader:
        read_name = row["read_name"]
        if read_name != prev_read_name:
            if not prev_nonrelevant_chim:
                count_virus_chims += 1
                for read in reads_collected:
                    writer.writerow(read)
            # reinit
            reads_collected = list()
            prev_nonrelevant_chim = False

            prev_read_name = read_name

        chr_donorA_is_virus = row["chr_donorA"] in viral_db_entries
        chr_donorB_is_virus = row["chr_acceptorB"] in viral_db_entries

        if (chr_donorA_is_virus ^ chr_donorB_is_virus) and "chrM" not in (row["chr_donorA"], row["chr_acceptorB"]):
            reads_collected.append(row)

        else:
            # ignore those that show up as human/human or virus/virus chimeric reads - not to be trusted as human/virus evidence
            prev_nonrelevant_chim = True

    # get last one
    if not prev_nonrelevant_chim:
        count_virus_chims += 1
        for read in reads_collected:
            writer.writerow(read)

    logger.info(f"Extracted {count_virus_chims} human/virus chimeric reads")

    fh.close()
    ofh.close()

    sys.exit(0)


if __name__ == "__main__":
    main()
