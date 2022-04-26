#!/usr/bin/env python

import pandas as pd
import sys

infile = sys.argv[-2]
output = sys.argv[-1]

ra = pd.read_table(infile, sep = "\t", header = 0)
wf = open(output + ".txt", "w")

i = 0
for row in ra.itertuples():
    i += 1
    r = ">r" + str(i)
    chr = row.chromosome
    pos = row.position
    cp = chr+":"+str(pos)
    ref = row.ref
    alt1 = row.alt1
    alt2 = row.alt2
    allele = sorted(set([ref, alt1, alt2]))
    if len(allele) > 2:
        allele1 = "".join([ref, alt1])
        wf.write(r + " " + cp + " " + allele1 + "\n")
        allele2 = "".join([ref, alt2])
        wf.write(r + " " + cp + " " + allele2 + "\n")
    else:
        allele = "".join(allele)
        wf.write(r + " " + cp + " " + allele + "\n")

wf.close()