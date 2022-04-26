#!/usr/bin/env python

"""
Merge all read count data from htseq-count output
"""

import pandas as pd
import os
import argparse

def args_parser():
    '''parser the arguments from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="Merge read count data from each sample and build a read count matrix. \n\
        Usage: htseq_merge.py -countdir <directory> -O <output>")
    parser.add_argument("-countdir", "--countdir", help = "the directory where has the read count files to be merged. ")
    parser.add_argument("-O", "--output", help = "assign the output name. ")
    args = parser.parse_args()
    return args

def merge_count(args):
    '''Require count data of ref or alt. '''
    cntdir = args.countdir
    out = args.output
    cnt_list = []
    htseqfiles = sorted([infile for infile in os.listdir(cntdir) if infile.endswith(".txt")])
    ### walk through all the files
    for infile in htseqfiles:
        pre = infile.split(".txt")[0]
        htseqcnt = pd.read_table(os.path.join(cntdir, infile), sep = '\t', header = None, comment = '_')
        htseqcnt = htseqcnt.rename(columns = {0:'gene', 1:pre})
        # htseqcnt['gene'] = htseqcnt['gene'].map(lambda x:x.lstrip("mRNA:"))
        # set gene column as index row 
        htseqcnt = htseqcnt.set_index("gene")
        cnt_list.append(htseqcnt)
    htseq = pd.concat(cnt_list, axis = 1)
    htseq.to_csv(out+ ".txt", sep = '\t', index = True)

def main():
    args = args_parser()
    merge_count(args)

##################
##### Run it #####
##################

if __name__ == "__main__":
    main()
