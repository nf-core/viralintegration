#!/usr/bin/env python
# encoding: utf-8

#########################################
# Author: [brianjohnhaas](https://github.com/brianjohnhaas)
# File: find_closest.py
# Source: https://github.com/broadinstitute/CTAT-VirusIntegrationFinder/blob/master/util/find_closest.py
# Source+commit: https://github.com/broadinstitute/CTAT-VirusIntegrationFinder/blob/677ddfc002009400ff5a0879ad2cf670933fdc66/util/find_closest.py
# Download Date: 2022-12-28, commit: 677ddfc
# This source code is licensed under the BSD 3-Clause license
#########################################

import sys
import argparse
import os
from subprocess import check_call

import pandas as pd

arguments = argparse.ArgumentParser(description="Add closest upstream and downstream genes")

arguments.add_argument(
    "-i",
    required=True,
    type=str,
    help="virus insertion finder summary tsv file",
)

arguments.add_argument(
    "-o",
    required=True,
    type=str,
    help="output tsv file",
)

arguments.add_argument(
    "--gtf",
    required=True,
    type=str,
    help="Host GTF file",
)

args = arguments.parse_args()

gtf = args.gtf
output_file = args.o
sorted_gtf = os.path.abspath("out.sorted.gtf")
input_tsv = args.i


input_tsv = pd.read_csv(input_tsv, sep="\t")
if len(input_tsv) == 0:
    print("no insertions reported.", file=sys.stderr)
    sys.exit(0)


with open(sorted_gtf, "wt") as f:
    check_call(["bedtools", "sort", "-i", gtf], stdout=f)

unique_chromosomes = set()
with open("out.sorted.gtf", "rt") as f:
    for line in f:
        if not line.startswith("#"):
            tokens = line.split("\t")
            unique_chromosomes.add(tokens[0].lower())

# convert tsv to bed
# entry	chrA	coordA	orientA	chrB	coordB	orientB	primary_brkpt_type	num_primary_reads	num_supp_reads	total

unsorted_bed_file = os.path.abspath("query.bed")
query_bed_file = os.path.abspath("query_sorted.bed")
upstream_gtf = os.path.abspath("upstream.gtf")
downstream_gtf = os.path.abspath("downstream.gtf")
with open(unsorted_bed_file, "wt") as f:
    for i in range(len(input_tsv)):
        entry = input_tsv.iloc[i]
        chra = str(entry["chrA"])
        chrb = str(entry["chrB"])
        if chra.lower() in unique_chromosomes:
            f.write(str(chra) + "\t" + str(entry["coordA"]) + "\t" + str(entry["coordA"]) + "\t" + str(i) + "\n")
        if chrb.lower() in unique_chromosomes:
            f.write(chrb + "\t" + str(entry["coordB"]) + "\t" + str(entry["coordB"]) + "\t" + str(i) + "\n")
with open(query_bed_file, "wt") as f:
    check_call(["bedtools", "sort", "-i", unsorted_bed_file], stdout=f)
with open(upstream_gtf, "wt") as f:
    # Ignore features in B that are downstream of features in A
    check_call(["bedtools", "closest", "-id", "-D", "a", "-a", query_bed_file, "-b", sorted_gtf], stdout=f)

with open(downstream_gtf, "wt") as f:
    # Ignore features in B that are downstream of features in A
    check_call(["bedtools", "closest", "-iu", "-D", "a", "-a", query_bed_file, "-b", sorted_gtf], stdout=f)


def parse_bedtools_output(path, column_name):
    closest_df = pd.read_csv(path, sep="\t", header=None)
    gene_names = []
    for i in range(len(closest_df)):
        entries = closest_df.iloc[i][12].split(";")
        gene_id = None
        gene_name = None
        for entry in entries:
            tokens = entry.strip().split(" ")
            if tokens[0] == "gene_id":
                value = tokens[1]
                gene_id = value[1 : len(value) - 1]
            elif tokens[0] == "gene_name":
                value = tokens[1]
                gene_name = value[1 : len(value) - 1]
        gene_names.append(gene_name if gene_name is not None else gene_id)
    closest_df["gene_name"] = gene_names
    closest_df = closest_df.drop_duplicates(subset=3)
    closest_df = closest_df.set_index(3)
    closest_df[column_name] = (
        closest_df["gene_name"].astype(str) + closest_df[10].astype(str) + " " + closest_df[13].astype(str)
    )
    return closest_df[[column_name]]


# closest_df[3] contains id from BED file
# closest_df[10] strand
# closest_df[12] contains gene_id "ENSG00000223972.5"; gene_type "transc...
# closest_df[13] contains distance
input_tsv = input_tsv.join(parse_bedtools_output(upstream_gtf, "upstream"))
input_tsv = input_tsv.join(parse_bedtools_output(downstream_gtf, "downstream"))
input_tsv.to_csv(output_file, sep="\t", index=False)
