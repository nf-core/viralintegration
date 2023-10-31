#!/usr/bin/env python3

#########################################
# Author: [brianjohnhaas](https://github.com/brianjohnhaas)
# File: chimJ_to_virus_insertion_candidate_sites.py
# Source: https://github.com/broadinstitute/CTAT-VirusIntegrationFinder/blob/master/util/chimJ_to_virus_insertion_candidate_sites.py
# Source+commit: https://github.com/broadinstitute/CTAT-VirusIntegrationFinder/blob/034594f47ccf63237eef3d305d07496fe6aeba44/util/chimJ_to_virus_insertion_candidate_sites.py
# Download Date: 2022-12-28, commit: 034594f
# This source code is licensed under the BSD 3-Clause license
#########################################

import sys, os, re
import argparse
import subprocess
import logging
from collections import defaultdict
import pandas as pd

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
        description="Defines candidate virus insertion sites based on chimeric reads",
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
        "--human_genome_aggregation_dist",
        type=int,
        required=False,
        default=500,
        help="distance around top chimeric event breakpoint for aggregating supporting reads on human genome",
    )

    parser.add_argument(
        "--viral_genome_aggregation_dist",
        type=int,
        required=False,
        default=500,
        help="distance around top chimeric event breakpoint for aggregating supporting reads on viral genome",
    )

    parser.add_argument(
        "--output_prefix",
        "-o",
        dest="output_prefix",
        type=str,
        required=True,
        help="output prefix",
    )

    parser.add_argument(
        "--remove_duplicates_flag",
        action="store_true",
        default=False,
        help="exclude duplicate alignments",
    )

    parser.add_argument(
        "--max_multi_read_alignments",
        type=int,
        default=1,  # unique chimeric read alignments
        help="max number of multimaps for chimeric read alignments to be considered as evidence (1=unique only)",
    )

    parser.add_argument("--debug", action="store_true", help="debug mode")

    ###########################
    # Parse input arguments
    ###########################
    args_parsed = parser.parse_args()
    # Create constants
    chimJ_filename = args_parsed.chimJ
    viral_db_fasta_filename = args_parsed.viral_db_fasta
    human_aggregation_dist = args_parsed.human_genome_aggregation_dist
    viral_aggregation_dist = args_parsed.viral_genome_aggregation_dist
    output_prefix = args_parsed.output_prefix
    max_multi_read_alignments = args_parsed.max_multi_read_alignments

    if args_parsed.debug:
        logger.setLevel(logging.DEBUG)

    remove_duplicates_flag = args_parsed.remove_duplicates_flag

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

    ################################
    # Get all Chrom-Virus pairings and store in 'genome_pair_to_evidence'
    # each line is a Chimeric_read object
    ################################
    # Junction type translations
    junction_type_encoding = {
        "-1": "Span",
        "0": "Split",
        "1": "Split",  # GT/AG",
        "2": "Split",
    }  # CT/AC" }
    # orientation conversions
    opposite_orientation = {"+": "-", "-": "+"}

    # dictionary to hold the pairings
    genome_pair_to_evidence = defaultdict(set)

    duplicate_read_catcher = set()

    # Read in the junction file
    df = pd.read_csv(chimJ_filename, sep="\t")
    df = df.astype({"read_name": "str"})  # in case they look like integers.

    # ~~~~~~~~~~~~~~~~~~~
    # multimapping reads
    # ~~~~~~~~~~~~~~~~~~~
    # - under defaults, we discard multimapping reads

    logger.info("## Mulitimapping reads")
    logger.info("-max multi read alignments allowed: {}".format(max_multi_read_alignments))

    hitcounts = df.groupby("read_name").size().to_dict()
    df["hitcount"] = df["read_name"].map(hitcounts)

    logger.info(f"-chim alignments before applying max multi setting: {df.shape[0]}")
    # Filter by multimapping
    df = df[df["hitcount"] <= max_multi_read_alignments]
    logger.info(f"-chim alignments AFTER applying max multi setting: {df.shape[0]}")

    # Convert the junction types
    df["junction_type"].replace(
        {
            -1: "Span",
            0: "Split",
            1: "Split",  # GT/AG",
            2: "Split",
        },
        inplace=True,
    )

    # ~~~~~~~~~~~~~~~~~~~~~~~~~
    # Check to see if chrA or chrB are in viral data base
    # ~~~~~~~~~~~~~~~~~~~~~~~~~
    logger.info("-filtering out human--human chimeric entries")
    df = df[(df["chr_donorA"].isin(viral_db_entries)) ^ (df["chr_acceptorB"].isin(viral_db_entries))]

    logger.info(f"Chimeric insertions involving viruses and human: {df.shape[0]}")

    # ~~~~~~~~~~~~~~~~~~~~~~~~~
    ## reorient so host (human) genome is always in the + reference orientation.
    # ~~~~~~~~~~~~~~~~~~~~~~~~~
    idx = (df["chr_donorA"].isin(viral_db_entries)) & (df["strand_acceptorB"] == "-") | (
        df["chr_acceptorB"].isin(viral_db_entries) & (df["strand_donorA"] == "-")
    )

    df.loc[
        idx,
        [
            "chr_donorA",
            "brkpt_donorA",
            "strand_donorA",
            "repeat_left_lenA",
            "start_alnA",
            "cigar_alnA",
            "chr_acceptorB",
            "brkpt_acceptorB",
            "strand_acceptorB",
            "repeat_right_lenB",
            "start_alnB",
            "cigar_alnB",
        ],
    ] = df.loc[
        idx,
        [
            "chr_acceptorB",
            "brkpt_acceptorB",
            "strand_acceptorB",
            "repeat_right_lenB",
            "start_alnB",
            "cigar_alnB",
            "chr_donorA",
            "brkpt_donorA",
            "strand_donorA",
            "repeat_left_lenA",
            "start_alnA",
            "cigar_alnA",
        ],
    ].values

    # adjust the orientation
    df.loc[idx, ["strand_donorA", "strand_acceptorB"]] = df.loc[idx, ["strand_donorA", "strand_acceptorB"]].replace(
        {"+": "-", "-": "+"}
    )

    if remove_duplicates_flag:
        logger.info("## Duplicate read alignments")
        # ~~~~~~~~~~~~~~~~~~~~~~~~~
        # Duplicates
        # ~~~~~~~~~~~~~~~~~~~~~~~~~
        # Add duplicated annotations
        # identify duplicated read names as well as duplicated insertions
        # duplicate insertions

        logger.info(f"Chimeric alignments before filter duplicates: {df.shape[0]}")

        df["Duplicate"] = df[
            [
                "chr_donorA",
                "start_alnA",
                "brkpt_donorA",
                "strand_donorA",
                "cigar_alnA",
                "chr_acceptorB",
                "start_alnB",
                "brkpt_acceptorB",
                "strand_acceptorB",
                "cigar_alnB",
            ]
        ].duplicated(keep="first")

        # df.to_csv("dups_marked.tsv", sep="\t")

        df = df[df["Duplicate"] == False]

        logger.info(f"Chimeric alignments AFTER filter duplicates: {df.shape[0]}")

    # Run through each line in the chimeric junctions file
    for idx, row in df.iterrows():
        # print(idx, row)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Create Dictionary holding read objects for each genome-virus pairing
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Output the chimeric read line as an object
        chim_read = Chimeric_read(
            row.chr_donorA,
            row.brkpt_donorA,
            row.strand_donorA,
            row.chr_acceptorB,
            row.brkpt_acceptorB,
            row.strand_acceptorB,
            row.junction_type,
            row.read_name,
        )
        # Creating Genome-Virus pairings
        # dictionary KEY
        genome_pair = "^".join(
            [
                row.chr_donorA,
                row.chr_acceptorB,
            ]
        )
        # Add the chimeric read lines to a dictionary with key as the Chrom-Virus grouping
        ## dictionary holds all chromosome virus groupings
        genome_pair_to_evidence[genome_pair].add(chim_read)

    #########################################
    # Gather all Chimeric events
    #########################################
    # empty list to hold all events at the end
    all_chim_events = list()

    # Loop over the different chimeric pairs (genome_pair_to_evidence dict keys)
    for genome_pair in genome_pair_to_evidence:
        # Get the Chimeric_read object for the given Chrom-Virus pairing
        chim_reads = genome_pair_to_evidence[genome_pair]

        logger.debug("Genome pair: " + genome_pair)

        chim_events = group_chim_reads_into_events(
            chim_reads, viral_db_entries, human_aggregation_dist, viral_aggregation_dist
        )

        # //TODO: can try to link up events into paired insertion points for single events.

        all_chim_events.extend(chim_events)

    ## prioritize by total read support.
    all_chim_events = sorted(all_chim_events, key=lambda x: (x.get_read_support()[2], str(x)), reverse=True)

    ## #####################
    ## generate final report.
    output_filename_full = output_prefix + ".full.tsv"
    output_filename_abridged = output_prefix + ".abridged.tsv"

    output_filename_detailed_abridged = output_prefix + ".abridged.detailed.tsv"
    ofh_detailed = open(output_filename_detailed_abridged, "wt")

    with open(output_filename_abridged, "wt") as ofh:
        with open(output_filename_full, "wt") as ofh_full:
            # print header
            ofh.write(
                "\t".join(
                    [
                        "entry",
                        "chrA",
                        "coordA",
                        "orientA",
                        "chrB",
                        "coordB",
                        "orientB",
                        "primary_brkpt_type",
                        "num_primary_reads",
                        "num_supp_reads",
                        "total",
                    ]
                )
                + "\n"
            )
            ofh_full.write(
                "\t".join(
                    [
                        "entry",
                        "chrA",
                        "coordA",
                        "orientA",
                        "chrB",
                        "coordB",
                        "orientB",
                        "primary_brkpt_type",
                        "num_primary_reads",
                        "num_supp_reads",
                        "total",
                        "readnames",
                    ]
                )
                + "\n"
            )

            ofh_detailed.write(
                "\t".join(
                    [
                        "parent_entry",
                        "entry",
                        "chrA",
                        "coordA",
                        "orientA",
                        "chrB",
                        "coordB",
                        "orientB",
                        "primary_brkpt_type",
                        "num_primary_reads",
                    ]
                )
                + "\n"
            )

            for chim_event in all_chim_events:
                chim_event.refine_insertion_coordinates()

                print(chim_event.get_event_accession() + "\t" + str(chim_event), file=ofh)
                supporting_reads = chim_event.get_readnames()
                print(
                    chim_event.get_event_accession()
                    + "\t"
                    + str(chim_event)
                    + "\t"
                    + ",".join(sorted(supporting_reads)),
                    file=ofh_full,
                )

                # detailed entry-level event report
                parent_chim_event = chim_event
                for event in [chim_event] + chim_event.chimeric_events_absorbed:
                    print(
                        parent_chim_event.get_event_accession()
                        + "\t"
                        + event.get_event_accession()
                        + "\t"
                        + event.get_coordstring()
                        + "\t"
                        + str(len(event.chimeric_reads_list)),
                        file=ofh_detailed,
                    )

    ofh_detailed.close()

    logger.info("-wrote output to {}".format(output_filename_abridged))

    sys.exit(0)


def group_chim_reads_into_events(chim_reads_list, viral_db_entries, human_aggregation_dist, viral_aggregation_dist):
    """
    Function to groups all chimeric events.
        if chimeric events occur in the same orientation and within a specific distance
        they will be combined

    Inputs;
    : chim_reads_list  : list of read objects
    : aggregation_dist : distance value
    Returns;
    : chim_events : a list of the chimertic events
    what are top event reads?
    """

    chim_events = list()

    remaining_reads = chim_reads_list

    # while reads still remain, run
    while remaining_reads:
        # Get the insertion event with the most insertions (counts)
        top_event_reads, remaining_reads = gather_top_event_reads(remaining_reads)
        chim_event = Chimeric_event(top_event_reads)

        # IF this chimeric event is the same orientation and within a specified distance from another event
        #   combine the event with the existing event in the chim_events list
        # IF NOT, add the new event to the chim_events list
        if not supplements_existing_event(
            chim_event, chim_events, viral_db_entries, human_aggregation_dist, viral_aggregation_dist
        ):
            # add new event to chim_events list
            chim_events.append(chim_event)
            logger.debug("-logging chimeric event: " + str(chim_event))

    return chim_events


def supplements_existing_event(
    chim_event, chim_events_list, viral_db_entries, human_aggregation_dist, viral_aggregation_dist
):
    """
    Function to check if the current chimeric read/event can be combined with an existing event
    """

    (agg_dist_A, agg_dist_B) = (
        (human_aggregation_dist, viral_aggregation_dist)
        if chim_event.chrB in viral_db_entries
        else (viral_aggregation_dist, human_aggregation_dist)
    )

    # Loop over stored Chimeric events
    for prev_chim_event in chim_events_list:
        # if the orientation is the same
        # if the read falls within the given distance from each other
        # then combine them

        if (
            chim_event.orientA == prev_chim_event.orientA
            and abs(chim_event.coordA - prev_chim_event.coordA) <= agg_dist_A
            and chim_event.orientB == prev_chim_event.orientB
            and abs(chim_event.coordB - prev_chim_event.coordB) <= agg_dist_B
        ):
            logger.debug("-adding {} as supplement to {}".format(str(chim_event), str(prev_chim_event)))
            # append this chimeric event to the previous
            prev_chim_event.absorb_nearby_chim_event(chim_event)

            return True

    return False


def gather_top_event_reads(reads_list):
    # count reads according to breakpoint.

    brkpt_counter = defaultdict(int)
    brkpt_type = dict()

    # Get break point types (split or span) into dictonary form
    #   example {'218269983^39325': 'Span'}
    # count the number of times this breakpoint is seen
    for read in reads_list:
        brkpt = "{}^{}".format(read.coordA, read.coordB)
        brkpt_counter[brkpt] += 1
        brkpt_type[brkpt] = read.splitType

    # prioritize split reads over spanning reads.
    priority = {"Split": 1, "Span": 0}

    # sort the breakpoint values
    # example output ['215417625^39519', '75063812^39394]
    # Sort forst by read type (split,span) then by read count
    #   Puts split first
    sorted_brkpts = sorted(
        brkpt_counter.keys(),
        key=lambda x: (priority[brkpt_type[x]], brkpt_counter[x], str(x)),
        reverse=True,
    )

    # Take the first value
    top_brkpt = sorted_brkpts[0]

    top_event_reads = list()
    remaining_reads = list()

    # Loop over reads to get the top breakpoint
    for read in reads_list:
        brkpt = "{}^{}".format(read.coordA, read.coordB)
        if brkpt == top_brkpt:
            top_event_reads.append(read)
        else:
            remaining_reads.append(read)

    return top_event_reads, remaining_reads


class Chimeric_read:
    """
    Object that parses the read
    """

    def __init__(self, chrA, coordA, orientA, chrB, coordB, orientB, splitType, readname=None):
        self.chrA = chrA
        self.coordA = coordA
        self.orientA = orientA
        self.chrB = chrB
        self.coordB = coordB
        self.orientB = orientB
        self.splitType = splitType
        self.readname = readname
        self._refined = False  # can only refine coodinates once!

    def __repr__(self):
        return "\t".join(
            [
                self.chrA,
                str(self.coordA),
                self.orientA,
                self.chrB,
                str(self.coordB),
                self.orientB,
                self.splitType,
            ]
        )


class Chimeric_event(Chimeric_read):
    """
    Object
    """

    def __init__(self, chimeric_reads_list):
        self.chimeric_reads_list = chimeric_reads_list

        self.chimeric_events_absorbed = (
            list()
        )  # for neighboring chim events that also support this primary event but are either spanning or split w/ different nearby brkpt

        example_read = chimeric_reads_list[0]
        super().__init__(
            example_read.chrA,
            example_read.coordA,
            example_read.orientA,
            example_read.chrB,
            example_read.coordB,
            example_read.orientB,
            example_read.splitType,
        )

    def absorb_nearby_chim_event(self, chim_event):
        self.chimeric_events_absorbed.append(chim_event)

    def get_read_support(self):
        num_chimeric_reads = len(self.chimeric_reads_list)

        num_absorbed_reads = 0
        for chim_event in self.chimeric_events_absorbed:
            num_absorbed_reads += len(chim_event.chimeric_reads_list)

        num_total_reads = num_chimeric_reads + num_absorbed_reads

        return num_chimeric_reads, num_absorbed_reads, num_total_reads

    def refine_insertion_coordinates(self):
        if self._refined:
            raise RuntimeError("chimeric event already refined, can only do this once")
        self._refined = True

        if self.splitType != "Span":
            # Split read, coodinates fixed.
            return

        # for spanning reads, use closest boundary based on orientation.
        coordA_vals = list()
        coordB_vals = list()
        coordA_vals.append(self.coordA)
        coordB_vals.append(self.coordB)

        for chim_event in self.chimeric_events_absorbed:
            coordA_vals.append(chim_event.coordA)
            coordB_vals.append(chim_event.coordB)

        self.coordA = max(coordA_vals) if self.orientA == "+" else min(coordA_vals)
        self.coordB = min(coordB_vals) if self.orientB == "+" else max(coordB_vals)

        return

    def __repr__(self):
        (
            num_chimeric_reads,
            num_absorbed_reads,
            num_total_reads,
        ) = self.get_read_support()

        return super().__repr__() + "\t{}\t{}\t{}".format(num_chimeric_reads, num_absorbed_reads, num_total_reads)

    def get_coordstring(self):
        return super().__repr__()

    def get_event_accession(self):
        return "~".join(
            [
                str(x)
                for x in (
                    self.chrA,
                    self.coordA,
                    self.orientA,
                    self.chrB,
                    self.coordB,
                    self.orientB,
                )
            ]
        )

    def get_readnames(self):
        readnames = list()

        for chim_event in [self] + self.chimeric_events_absorbed:
            for chim_read in chim_event.chimeric_reads_list:
                readnames.append(chim_read.readname)

        return readnames


if __name__ == "__main__":
    main()
