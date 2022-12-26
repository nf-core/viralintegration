#!/usr/bin/env python3

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
