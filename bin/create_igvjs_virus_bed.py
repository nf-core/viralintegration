#!/usr/bin/env python3
# encoding: utf-8
if __name__ == "__main__":
    import argparse
    import pandas as pd
    import sys

    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", type=str, required=True, help="prefix.virus_read_counts_summary.tsv")
    parser.add_argument("--output_prefix", type=str, required=True, help="Prefix for output files")
    parser.add_argument("--num_top_viruses", type=int, required=False, default=None, help="num top viruses")
    parser.add_argument(
        "--min_bases_covered", type=int, default=200, help="minimum number of viral genomic bases covered"
    )
    args = parser.parse_args()

    df = pd.read_csv(args.summary, sep="\t")  # virus	seqlen	mapped	chim_reads
    df = df[df["mapped"] > 0]
    df = df[df["n_bases_covered"] >= args.min_bases_covered]

    # df.sort_values(by='mapped', ascending=False, inplace=True) # no, keep incoming sort order.

    if args.num_top_viruses is not None and args.num_top_viruses >= 1 and len(df) > args.num_top_viruses:

        df = df.head(args.num_top_viruses)

    df.to_csv(args.output_prefix + ".igvjs.table.tsv", index=False, sep="\t")

    # bed file

    df["start"] = 0
    df.to_csv(
        args.output_prefix + ".igvjs.bed", index=False, header=False, columns=["virus", "start", "seqlen"], sep="\t"
    )

    sys.exit(0)
