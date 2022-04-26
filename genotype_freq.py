#!/usr/bin/env python

import pandas as pd
import argparse
import os
import numpy as np
from mpi4py import MPI

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="Filter bad SNPs by combining all samples' read-count data. \n \
        Usage: mpiexec -n <cores> python genotype_caller.py -dir <dir> -O <prefix>")
    parser.add_argument("-dir", "--directory", help = "the directory which contains all the read-count files. ")
    parser.add_argument("-O", "--output", help = "output file with bad SNPs  ")
    args = parser.parse_args()
    return args

def site_genotype(args, worker_tasks, rank):
    dir = args.directory
    i = 0
    file_list = []
    for f in worker_tasks[rank]:
        core_files = len(worker_tasks[rank])
        pre = f.split(".txt")[0]
        i += 1
        print(f, "Read the " + str(i) + " file in total " + str(core_files) + " on the core " + str(rank))
        count_df = pd.read_table(os.path.join(dir, f), header = 0, sep = '\t')
        ### informative SNP
        count_df = count_df[((count_df["alt1_count"] > 10) | (count_df["alt2_count"] > 10)) & (count_df['others']/(count_df['alt1_count'] + count_df['alt2_count'] + count_df['others'])*100 < 1)]
        ### assign genotype
        count_df['alt2_perc'] = count_df['alt2_count']/(count_df['alt1_count'] + count_df['alt2_count'])*100
        count_df[pre] = np.where((count_df['alt2_count'] < 10) & (count_df['alt2_perc'] < 5), "homo", "hete")
        # combine all genotyped sites
        count_gen = count_df[['chromosome', "position", pre]]
        count_gen = count_gen.set_index(keys = ["chromosome", "position"])
        file_list.append(count_gen)
    return file_list

def main():
    ### parallele run
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    ### parse arguments
    args = args_parser()
    output = args.output
    count_files = sorted([x for x in os.listdir(args.directory) if x.endswith(".txt")])
    if rank == 0:
        dir_num = len(count_files)
        worker_tasks = {w:[] for w in range(size)}
        w_idx = 0
        for i in range(dir_num):
            worker_tasks[w_idx].append(count_files[i])
            w_idx = (w_idx + 1) % size
        print("split work successfully!") # rank:files dictionary
    else:
        worker_tasks = None
    worker_tasks = comm.bcast(worker_tasks, root = 0)
    ### assign genotype for each site in each file
    site_files = site_genotype(args, worker_tasks, rank) # file list
    df = pd.concat(site_files, axis = 1) # concatenate for columns
    df.to_csv("rank_"+str(rank)+".txt", sep = "\t", index = True)
    df = comm.gather(df, root = 0) # return list
    if rank == 0:
        print("Combine all results on all cores!")
        df = pd.concat(df, axis = 1)
        df['homo'] = (df == "homo").sum(axis = 1)
        df['hete'] = (df == "hete").sum(axis = 1)
        df_g = df[['homo', 'hete']]
        ###
        df_g.to_csv(output+".txt", sep = "\t", index = True)
        print("Collect genotype frequency for each SNP site!")
        if output + ".txt" in os.listdir():
            for s in range(size):
                os.remove("rank_"+str(s)+".txt")
                print("remove all temporary files. ")

##############
### Run it ###
##############

if __name__ == "__main__":
    main()
