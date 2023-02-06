#!/usr/bin/env python3

#########################################
# Author: [brianjohnhaas](https://github.com/brianjohnhaas)
# File: region_gtf_to_bed.py
# Source: https://github.com/broadinstitute/CTAT-VirusIntegrationFinder/blob/master/util/region_gtf_to_bed.py
# Source+commit: https://github.com/broadinstitute/CTAT-VirusIntegrationFinder/blob/01be10e1041301cd5df4350664380d4563df1704/util/region_gtf_to_bed.py
# Download Date: 2022-12-28, commit: 01be10e
# This source code is licensed under the BSD 3-Clause license
#########################################

import sys, os, re


def main():
    if len(sys.argv) < 2:
        print(
            "\n\n\tusage: {} regions.gtf >  regions.bed\n\n".format(sys.argv[0]),
            file=sys.stderr,
        )
        sys.exit(1)

    regions_gtf_filename = sys.argv[1]

    with open(regions_gtf_filename, "r") as fh:
        for line in fh:
            line = line.rstrip()
            vals = line.split("\t")
            if len(vals) < 9:
                continue

            contig = vals[0]
            lend = vals[3]
            rend = vals[4]
            orient = vals[6]
            annot = vals[8]

            print("\t".join([contig, str(int(lend) - 1), rend, annot, ".", orient]))

    sys.exit(0)


if __name__ == "__main__":
    main()
