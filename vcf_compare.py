#!/usr/bin/env python

import pysam
import argparse
from Bio import SeqIO
import subprocess
import os 
import pandas as pd

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="Filter variants and rewrite the passed the variants into a new vcf file. \n\
    Usage: vcf_pass.py -vcf <snp.vcf> -R <reference.fa> -O <vcf_output>")
    parser.add_argument("-vcf1", "--vcf1", help="the first VCF file with PASS flag. ")
    parser.add_argument("-vcf2", "--vcf2", help = "the second VCF file with PASS flag. ")
    parser.add_argument("-R", "--reference", help="reference genome in fasta file")
    parser.add_argument("-O", "--output", help="output prefix")
    args=parser.parse_args()
    return args

def fasta_dict(args):
    '''parser reference genome and build dictionary'''
    ref_fasta = args.reference
    seq_parse=SeqIO.parse(ref_fasta, "fasta")
    seq_dict=SeqIO.to_dict(seq_parse)
    seq_dict_mut={}
    for entry in seq_dict:
        seq_dict_mut[seq_dict[entry].id] = str(seq_dict[entry].seq)
    return seq_dict_mut

def snp_collect(args, seq_dict):
    # pysam.VariantFile read vcf file and write header into new file
    output = args.output
    #df1 = pd.DataFrame(columns = ["chromosome", "position", "ref", "alt1"])
    df1 = open(output + "_1.txt", "w")
    df1.write("chromosome\tposition\tref\talt\n")
    df2 = open(output + "_2.txt", "w")
    df2.write("chromosome\tposition\tref\talt\n")
    # chromosome/contig with snp
    vcf1 = args.vcf1
    vcf2 = args.vcf2
    chr1 = {v.chrom for v in pysam.VariantFile(vcf1).fetch()}
    chr2 = {v.chrom for v in pysam.VariantFile(vcf2).fetch()}
    ###
    for c in sorted(chr1):
        seq_length = len(seq_dict[c])
        for v in pysam.TabixFile(vcf1, parser = pysam.asVCF()).fetch(c, 0, seq_length):
            pos = v.pos + 1
            alt = v.alt
            ref = v.ref
            df1.write("{0:s}\t{1:d}\t{2:s}\t{3:s}\n".format(c, pos, ref, alt))
    for c in sorted(chr2):
        seq_length = len(seq_dict[c])
        for v in pysam.TabixFile(vcf2, parser = pysam.asVCF()).fetch(c, 0, seq_length):
            pos = v.pos + 1
            alt = v.alt
            ref = v.ref
            df2.write("{0:s}\t{1:d}\t{2:s}\t{3:s}\n".format(c, pos, ref, alt))
    df1.close()
    df2.close()

def merge_snp(args):
    output = args.output
    df1 = pd.read_table(output + "_1.txt", sep = "\t", header = 0)
    df2 = pd.read_table(output + "_2.txt", sep = "\t", header = 0)
    df = df1.merge(df2, left_on = ["chromosome","position","ref"], right_on = ["chromosome","position","ref"], suffixes = ("_1", "_2"), how = "outer")
    df = df.sort_values(by = ["chromosome", "position"])
    #df.to_csv(output + ".txt", sep = "\t", index = False)
    df["ALT1"] = df["alt_1"].fillna(df["ref"])
    df["ALT2"] = df["alt_2"].fillna(df["ref"])
    df = df[["chromosome", "position", "ref", "ALT1", "ALT2"]]
    #df.to_csv(output + ".merge.txt", sep = "\t", index = False)
    df = df[df["ALT1"] != df["ALT2"]]
    df.to_csv(output + ".txt", sep = "\t", index = False)

def main():
    args = args_parser()
    seq_dict = fasta_dict(args)
    snp_collect(args, seq_dict)
    merge_snp(args)

############################
######### Run it ###########
############################

if __name__ == "__main__":
    main()
