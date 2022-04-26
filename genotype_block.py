#!/usr/bin/env python

import pandas as pd
import argparse
import numpy as np
import os

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description=" \n \
        Usage: genotype_block.py -chr <chr> -C <cnt.txt>")
    parser.add_argument("-chr", '--chromosome', help = 'The interested chromosomes you want to check, provide in a file. ')
    parser.add_argument("-C", "--count", help = "Read count file for the genotyping. ")
    parser.add_argument("-chrLen", '--chrLen', help = "chromosome length ")
    parser.add_argument("-O", "--output", default = 'default', help = "Prefix for the output file. ")
    args = parser.parse_args()
    return args

def outprefix(args):
    out = args.output
    c = args.count
    c = os.path.basename(c)
    pre = c.split(".txt")[0]
    if out == 'default':
        out = pre
    return out

def input_chr(args):
    ### genotype chromosome for any of them with informative SNPs
    chr = pd.read_table(args.chromosome, sep = "\t", header = 0)
    chr1 = set(chr.chromosome)
    cnt = pd.read_table(args.count, sep = '\t', header = 0)
    cnt = cnt[((cnt['alt1_count'] > 10) | (cnt['alt2_count'] > 10)) & (cnt['others']/(cnt['alt1_count'] + cnt['alt2_count'] + cnt['others'])*100 < 1)]
    chr2 = set(cnt.chromosome)
    chr = chr1.intersection(chr2)
    return chr

def genotyping(args, out, chr):
    chrlen = pd.read_table(args.chrLen, sep = "\t", header = 0)
    cnt = pd.read_table(args.count, sep = '\t', header = 0)
    cnt = cnt.sort_values(by = ["chromosome", "position"])
    cnt = cnt[((cnt['alt1_count'] > 10) | (cnt['alt2_count'] > 10)) & (cnt['others']/(cnt['alt1_count'] + cnt['alt2_count'] + cnt['others'])*100 < 1)]
    cnt['alt2_perc'] = cnt['alt2_count']/(cnt['alt1_count'] + cnt['alt2_count'])*100
    cnt['genotype'] = np.where((cnt['alt2_perc'] < 5) & (cnt['alt2_count'] < 10), 'homo', 'hete')
    cnt.to_csv(out + ".cnt", sep = '\t', index = False) # for further check
    break_record = open(out + ".brk", "w")
    break_record.write("chromosome\tgenotype\tstart\tend\ttransition_site\tlength\tnoise\n")
    ###
    for r in sorted(chr):
        #print(r)
        rlen = chrlen.loc[(chrlen['chromosome'] == r), 'length'].item()
        noise = []
        t1 = ''
        pos1 = ''
        sr = 1
        subcnt = cnt[cnt['chromosome'] == r]
        snp = len(subcnt)
        start_pos = min(subcnt.position)
        end_pos = max(subcnt.position)
        ti = subcnt.loc[(subcnt['position'] == start_pos), 'genotype'].item()
        ###
        nis = 0
        tsnp = 0
        for row in subcnt.itertuples():
            tsnp += 1
            if len(t1) > 0:
                ti = t1
                posi = pos1
            t1 = row.genotype
            pos1 = row.position # information for current row
            if (t1 != ti): # potential break
                #print('break', posi, pos1)
                if break_bh(subcnt, pos1, t1):
                    bp = (posi + pos1)//2
                    print("({}-{}) transition on {}".format(posi, pos1, r))
                    l = bp - sr + 1
                    pc = round((nis/tsnp) * 100, 2)
                    break_record.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(r, ti, sr, bp, str(posi) + "-" + str(pos1), l, str(pc)+ "%(" + str(nis) +"/" + str(tsnp) +")"))
                    ###
                    sr = bp + 1
                    nis = 0
                    tsnp = 0
                else:
                    nis += 1
                    subcnt.loc[(subcnt.position == pos1), 'genotype'] = ti
                    t1 = ti
                    noise.append(pos1)
            if row.position == end_pos: # capture the last row
                l = rlen - sr + 1
                pc = round((nis/tsnp) * 100, 2)
                break_record.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(r, t1, sr, rlen, '-', l, str(pc)+ "%(" + str(nis) +"/" + str(tsnp) +")"))
                #print("From {} to the end of {}".format(sr, r))
    break_record.close()

def break_bh(subcnt, pos1, t1):
    break_check = False
    checnt = subcnt[subcnt['position'] >= pos1]
    ###
    if (len(checnt) > 20) & ((max(checnt.position)-min(checnt.position)) > 50000):
        i = 1
        while True:
            #print(checnt)
            geno_cnt = subcnt[(subcnt['position'] >= pos1) & (subcnt['position'] < (pos1 + i*1000))]
            #print(geno_cnt)
            i += 1
            re = len(geno_cnt)
            rng = max(geno_cnt.position) - min(geno_cnt.position)
            #print(subcnt[(subcnt['position'] >= pos1) & (subcnt['position'] < (pos1 + 5000))])
            if (re > 20) and (rng > 50000):
                break
        #print("i", i, geno_cnt)
        ###
        #print(geno_cnt)
        t_num = list(geno_cnt.genotype).count(t1)
        t_perc = (t_num/len(list(geno_cnt.genotype))) * 100 # the percentage of 
        #print(t_perc, '%, line_num', t_num)
        if (t_perc >= 80): # based on the variants percentage around 0.8%
            break_check = True
    return break_check

def main():
    args = args_parser()
    out = outprefix(args)
    print(out)
    chr = input_chr(args)
    genotyping(args, out, chr) # fixed homozygous to heterozygous junction


##################
##### Run it #####
##################

if __name__ == "__main__":
    main()