#!/usr/bin/env python

import sys, os, re
import pandas as pd
import subprocess as sp

if len(sys.argv) < 4:
    exit("usage: {} virus_read_counts.tsv input_virus.fa  output_virus.fa\n\n".format(sys.argv[0]))

virus_igvjs_bed = sys.argv[1]
input_virus_fa = sys.argv[2]
output_virus_fa = sys.argv[3]

if os.path.exists(output_virus_fa):
    exit("Error, output file {} already exists.  Remove or relocate it before proceeding".format(output_virus_fa))


sp.check_call("touch {}".format(output_virus_fa), shell=True)  # just in case, avoid missing outfile if no virus

df = pd.read_csv(virus_igvjs_bed, sep="\t", names=["virus", "start", "end"])
for _, row in df.iterrows():
    virus_acc = row["virus"]
    cmd = 'samtools faidx {} "{}" >> {}'.format(input_virus_fa, virus_acc, output_virus_fa)
    print(cmd, file=sys.stderr)
    sp.check_call(cmd, shell=True)

sys.exit(0)
