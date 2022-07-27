#!/usr/bin/env python3
# encoding: utf-8


import argparse
import csv
import json

import sys


arguments = argparse.ArgumentParser(
    prog="Virus Insertion Site viewer JSON Maker",
    description="Makes a JSON file for a directory of results from Virus Insertion Finder",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)


arguments.add_argument(
    "--VIF_summary_tsv",
    required=True,
    type=str,
    help="virus insertion finder summary tsv file",
)

arguments.add_argument(
    "--json_outfile", required=True, type=str, help="The output json file to create",
)


args = arguments.parse_args()

vif_summary_tsv_filename = args.VIF_summary_tsv
json_outfile = args.json_outfile


# Dict to be translated to JSON object
dict_json = {}
dict_json["fusions"] = []

vif_info_tokens = [
        "contig",
        "chrA",
        "coordA",
        "orientA",
        "chrB",
        "coordB",
        "orientB",
        "split",
        "span",
        "total",
        "upstream",
        "downstream"
]

# Make fusion detail
with open(vif_summary_tsv_filename, "r") as fh:
    csv_reader = csv.DictReader(fh, delimiter="\t")

    for row in csv_reader:
        dict_json["fusions"].append(row)

# Store as a json object
with open(json_outfile, "w") as write_json:
    write_json.write(json.dumps(dict_json, sort_keys=True, indent=2))

sys.exit(0)
