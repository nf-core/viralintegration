#!/usr/bin/env python

#########################################
# Author: [brianjohnhaas](https://github.com/brianjohnhaas)
# File: create_igvjs_virus_fa.py
# Source: https://github.com/broadinstitute/CTAT-VirusIntegrationFinder/blob/master/util/create_igvjs_virus_fa.py
# Source+commit: https://github.com/broadinstitute/CTAT-VirusIntegrationFinder/blob/e7fed4354c3eee1fac3accaac795c0dab5f9afd1/util/create_igvjs_virus_fa.py
# Download Date: 2022-12-28, commit: e7fed43
# This source code is licensed under the BSD 3-Clause license
#########################################

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
