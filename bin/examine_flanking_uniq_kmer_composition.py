#!/usr/bin/env python3

#########################################
# Author: [brianjohnhaas](https://github.com/brianjohnhaas)
# File: examine_flanking_uniq_kmer_composition.py
# Source: https://github.com/broadinstitute/CTAT-VirusIntegrationFinder/blob/master/util/examine_flanking_uniq_kmer_composition.py
# Source+commit: https://github.com/broadinstitute/CTAT-VirusIntegrationFinder/blob/5b6d71b3d0e4a83a4b41dfb0a0b7ab2c8b18ef83/util/examine_flanking_uniq_kmer_composition.py
# Download Date: 2022-12-28, commit: 5b6d71b
# This source code is licensed under the BSD 3-Clause license
#########################################

import sys, os, re
import pandas as pd
import argparse


def main():
    parser = argparse.ArgumentParser(
        description="include fraction unique kmer content metric",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("--vif_tsv", type=str, required=True, help="vif tsv input file")
    parser.add_argument("--output", type=str, required=True, help="vif output including the kmer metrics")
    parser.add_argument(
        "--min_frac_uniq",
        type=float,
        required=False,
        default=0.0,
        help="minimum fraction unique for flanking region kmer content k=5",
    )

    args = parser.parse_args()

    vif_tsv_filename = args.vif_tsv
    output_filename = args.output
    min_frac_uniq = args.min_frac_uniq

    df = pd.read_csv(vif_tsv_filename, sep="\t")

    df = df.apply(examine_unique_kmer_fraction, axis=1)  # fU = fraction unique

    df = df[(df["flankA_fU"] >= min_frac_uniq) & (df["flankB_fU"] >= min_frac_uniq)]

    df["flankA_fU"] = df["flankA_fU"].apply(lambda x: "{:.3f}".format(x))
    df["flankB_fU"] = df["flankB_fU"].apply(lambda x: "{:.3f}".format(x))

    df.to_csv(output_filename, sep="\t", index=False)

    sys.exit(0)


def examine_unique_kmer_fraction(row):
    row["flankA_fU"] = fraction_unique(row["flankA"])
    row["flankB_fU"] = fraction_unique(row["flankB"])

    return row


def fraction_unique(nuc_seq):
    K = 5
    uniq_kmers = set()
    kmer_count = 0
    nuc_seq = nuc_seq.upper()
    for i in range(0, len(nuc_seq) - K):
        kmer = nuc_seq[i : i + K]
        uniq_kmers.add(kmer)
        kmer_count += 1

    frac_uniq = len(uniq_kmers) / kmer_count

    return frac_uniq


if __name__ == "__main__":
    main()
