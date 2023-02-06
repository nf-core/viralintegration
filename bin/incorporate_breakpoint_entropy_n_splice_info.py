#!/usr/bin/env python

#########################################
# Author: [brianjohnhaas](https://github.com/brianjohnhaas)
# File: incorporate_breakpoint_entropy_n_splice_info.py
# Source: https://github.com/broadinstitute/CTAT-VirusIntegrationFinder/blob/master/util/incorporate_breakpoint_entropy_n_splice_info.py
# Source+commit: https://github.com/broadinstitute/CTAT-VirusIntegrationFinder/blob/6267c557034ba4a7aa5b22f5c20bf69eff85031c/util/incorporate_breakpoint_entropy_n_splice_info.py
# Download Date: 2022-12-28, commit: 6267c55
# This source code is licensed under the BSD 3-Clause license
#########################################

import sys, os, re
import pysam
from collections import defaultdict
import logging
import argparse
import pandas as pd
import subprocess
import math


logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description="add breakpoint entropy stats", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("--vif_tsv", required=True, type=str, help="VIF full tsv containing evidence read names")
    parser.add_argument("--ref_genome_fasta", required=True, type=str, help="ref genome fasta file")
    parser.add_argument("--viral_genome_fasta", required=True, type=str, help="viral genome fasta file")
    parser.add_argument("--output", required=True, type=str, help="output tsv file containing additional entropy stats")

    args = parser.parse_args()

    vif_tsv_filename = args.vif_tsv
    ref_genome_fasta_filename = args.ref_genome_fasta
    viral_genome_fasta_filename = args.viral_genome_fasta
    output_filename = args.output

    vif_df = pd.read_csv(vif_tsv_filename, sep="\t")
    print(vif_df.head())

    if len(vif_df) == 0:
        logger.info("-no candidates to pursue")
        subprocess.check_call(f"cp {vif_tsv_filename} {output_filename}", shell=True)
        sys.exit(0)

    ref_genome_fai_filename = ref_genome_fasta_filename + ".fai"
    if not os.path.exists(ref_genome_fai_filename):
        run_cmd("samtools faidx {}".format(ref_genome_fasta_filename))

    viral_genome_fai_filename = viral_genome_fasta_filename + ".fai"
    if not os.path.exists(viral_genome_fai_filename):
        run_cmd("samtools faidx {}".format(viral_genome_fasta_filename))

    virus_accs = set(pd.read_csv(viral_genome_fai_filename, sep="\t", header=None)[0].tolist())

    flank_len = 30

    def get_breakpoint_flanking_seq(acc, coord, orient, left_or_right_side):
        # coord represents donor or acceptor site position
        # so need to adjust for sequence extraction based on left/right and orientation info.

        if left_or_right_side == "left":
            if orient == "+":
                anchor_left = coord - 1 - flank_len + 1
                anchor_right = coord - 1 + 2
            elif orient == "-":
                anchor_left = coord + 1 - 2
                anchor_right = coord + 1 + flank_len - 1

            else:
                raise RuntimeError(f"cannot recognize orient {orient}")

        elif left_or_right_side == "right":
            if orient == "+":
                anchor_left = coord + 1 - 2
                anchor_right = coord + 1 + flank_len - 1
            elif orient == "-":
                anchor_left = coord - 1 - flank_len + 1
                anchor_right = coord - 1 + 2
            else:
                raise RuntimeError(f"cannot recognize orient {orient}")
        else:
            raise RuntimeError(f"cannot recognize left_or_right_side {left_or_right_side}")

        if acc in virus_accs:
            seq_region = extract_seqrange(acc, viral_genome_fasta_filename, anchor_left, anchor_right)
        else:
            seq_region = extract_seqrange(acc, ref_genome_fasta_filename, anchor_left, anchor_right)

        if orient == "-":
            seq_region = revcomp(seq_region)

        seq_region = list(seq_region.lower())

        if left_or_right_side == "left":
            seq_region[-2] = seq_region[-2].upper()
            seq_region[-1] = seq_region[-1].upper()
        else:
            seq_region[0] = seq_region[0].upper()
            seq_region[1] = seq_region[1].upper()

        seq_region = "".join(seq_region)

        return seq_region

    vif_df["flankA"] = vif_df.apply(
        lambda row: get_breakpoint_flanking_seq(row["chrA"], row["coordA"], row["orientA"], "left"), axis=1
    )
    vif_df["flankB"] = vif_df.apply(
        lambda row: get_breakpoint_flanking_seq(row["chrB"], row["coordB"], row["orientB"], "right"), axis=1
    )

    vif_df["entropyA"] = vif_df["flankA"].apply(lambda x: compute_entropy(x, "left"))
    vif_df["entropyB"] = vif_df["flankB"].apply(lambda x: compute_entropy(x, "right"))

    # simplify formatting.
    vif_df["entropyA"] = vif_df["entropyA"].apply(lambda x: "{:.3f}".format(x))
    vif_df["entropyB"] = vif_df["entropyB"].apply(lambda x: "{:.3f}".format(x))

    vif_df["splice_type"] = vif_df.apply(lambda row: get_splice_info(row["flankA"], row["flankB"]), axis=1)

    logger.info("-writing outputfile: {}".format(output_filename))
    vif_df.to_csv(output_filename, sep="\t", index=False)

    logger.info("-done")

    sys.exit(0)


def extract_seqrange(acc, fasta_filename, lend, rend):
    cmd = f"samtools faidx {fasta_filename} {acc}:{lend}-{rend}"

    logger.info(cmd)

    seq_txt = subprocess.check_output(cmd, shell=True).decode()
    seq_txt = "".join(seq_txt.split("\n")[1:])

    return seq_txt


def compute_entropy(seq_txt, left_or_right):
    seq_txt = seq_txt.upper()

    # remove splice dinucs
    if left_or_right == "right":
        seq_txt = seq_txt[2:]
    elif left_or_right == "left":
        seq_txt = seq_txt[:-2]
    else:
        raise RuntimeError(f"Error, not recognizing left_or_right [{left_or_right}] as left or right")

    char_counter = defaultdict(int)
    for char in seq_txt:
        char_counter[char] += 1

    num_chars = len(seq_txt)
    entropy = 0.0
    for char, count in char_counter.items():
        p = count / num_chars
        entropy += p * math.log2(1 / p)

    return entropy


revcomp_translation = {"G": "C", "g": "c", "C": "G", "c": "g", "A": "T", "a": "t", "T": "A", "t": "a"}


def revcomp(sequence):
    sequence = list(sequence)
    sequence = sequence[::-1]  # rev it

    sequence = [revcomp_translation[x] if x in revcomp_translation else x for x in sequence]

    sequence = "".join(sequence)

    return sequence


def get_splice_info(flankA, flankB):
    splice_signals = {"GT-AG", "GC-AG"}

    spliceA = flankA[-2:]
    spliceB = flankB[0:2]

    splice_candidate = f"{spliceA}-{spliceB}"

    if splice_candidate in splice_signals:
        return splice_candidate
    else:
        rev_splice_candidate = revcomp(splice_candidate)
        if rev_splice_candidate in splice_signals:
            return rev_splice_candidate

    return "."  # just a placeholder to indicate no canonical splice signal identified.


def run_cmd(cmd):
    logger.info("CMD: {}".format(cmd))
    subprocess.check_call(cmd, shell=True)
    return


if __name__ == "__main__":
    main()
