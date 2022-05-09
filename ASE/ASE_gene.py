

import pandas as pd
import pysam
import argparse
from mpi4py import MPI
import os

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="Takes SNPs file and perform allele-specific read count on gene-basis. \n\
        Usage: mpiexec -n <cores> allele_count.py -SNP <SNP> -gtf <GTF> -bam <BAM> -O <prefix>")
    parser.add_argument("-SNP", "--SNP", help = "SNP on exonic region with position on 1-basis. ")
    parser.add_argument("-gtf", "--gtf", help="sorted gtf file with its tabix index file. ")
    parser.add_argument("-bam", "--bam", help = "the sorted bam with index file in the same directory. ")
    parser.add_argument("-O", "--output", help="output table with allele-specific read count on gene-basis. ")
    args = parser.parse_args()
    return args

def gene_pos(args):
    '''acquire genes position information in the gtf file. '''
    gene_dict = {}
    gtf_handle = pysam.TabixFile(args.gtf, parser = pysam.asGTF())
    contig_set = sorted({ginfo.contig for ginfo in gtf_handle.fetch()})
    for contig in contig_set:
        gene_dict[contig] = {gr.gene_id.split(":")[-1]:[gr.start, gr.end] for gr in gtf_handle.fetch(contig) if gr.feature == "gene"} # 0-based start position
    return gene_dict

def gene_strand(args, gene_dict):
    '''gene oritation based on gtf file. '''
    gstrand = {}
    gtf_handle = pysam.TabixFile(args.gtf, parser = pysam.asGTF())
    # all genes
    for contig in gene_dict:
        for gid in gene_dict[contig]:
            start = gene_dict[contig][gid][0] # 0-based 
            end = gene_dict[contig][gid][-1] # 1-based
            for ginfo in gtf_handle.fetch(contig, start, end):
                if (ginfo.gene_id.split(":")[-1] == gid and ginfo.feature == "gene"):
                    gstrand[gid] = ginfo.strand
    return gstrand

def read_strand(read, gstrand, gene_id):
    if gstrand[gene_id] == "+":
        if ((read.is_read1 and read.is_reverse) or (read.is_read2 and read.is_reverse == False)):
            return True
        else:
            return False
    elif gstrand[gene_id] == "-":
        if ((read.is_read1 and read.is_reverse == False) or (read.is_read2 and read.is_reverse)):
            return True
        else:
            return False
    else:
        return False

def variants_allele(args, gid, var_df, i, gene_dict, gstrand, bam_handle):
    '''allele-specific read count for each gene. '''
    ### container
    gene_read = {}
    alt1_read = {}
    alt2_read = {}
    non_read = {}
    ### ouput table for allele-specific read count
    # var_table = pd.read_table(args.SNP, sep = "\t", header = 0)
    vid = gid[i]
    sub_var = var_df[var_df["gene_id"] == vid]
    ctg = sorted(set(sub_var["chromosome"]))[0]
    ###
    gsta = gene_dict[ctg][vid][0] # 0-based
    gend = gene_dict[ctg][vid][-1] # 1-based
    #grange = "-".join([str(gsta+1), str(gend)])
    var_num = len(sub_var) # SNPs number on particular gene
    ### container for read labels
    gene_read[vid] = set()
    alt1_read[vid] = set()
    alt2_read[vid] = set()
    non_read[vid] = set()
    ### go through each SNP on particular gene
    for row in sub_var.itertuples():
        vpos = row.position # 1-based
        alt1 = row.alt1
        alt2 = row.alt2
        for read_column in bam_handle.pileup(ctg, vpos-1, vpos): # 0-based
            if read_column.reference_pos == vpos-1:
                for bam_read in read_column.pileups:
                    read_name = bam_read.alignment.query_name
                    if (("NH", 1) in bam_read.alignment.tags and bam_read.query_position != None and read_strand(bam_read.alignment, gstrand, vid)):
                        ###
                        gene_read[vid].add(read_name)
                        ###
                        allele = bam_read.alignment.query_sequence[bam_read.query_position]
                        if allele == alt1:
                            alt1_read[vid].add(read_name)
                        elif allele == alt2:
                            alt2_read[vid].add(read_name)
                        else:
                            non_read[vid].add(read_name)
    ### check if reads with allele in conflict
    conflict = alt1_read[vid].intersection(alt2_read[vid])
    cnft = len(conflict)
    for i in conflict:
        alt1_read[vid].remove(i)
        alt2_read[vid].remove(i)
    ###
    gene_read_count = len(gene_read[vid])
    alt1_read_count = len(alt1_read[vid])
    alt2_read_count = len(alt2_read[vid])
    non_read_count = len(non_read[vid])
    return vid, gene_read_count, alt1_read_count, alt2_read_count, non_read_count, cnft

def check_rm(args, size):
    output = args.output
    if output + ".txt" in os.listdir():
        for s in range(size):
            os.remove(output + "_tmp" + str(s) + ".txt")

def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    ###
    args = args_parser()
    gene_dict = gene_pos(args)
    gstrand = gene_strand(args, gene_dict)
    bam_handle = pysam.AlignmentFile(args.bam, "rb")
    ###
    if rank == 0:
        print("Build gene position dictionary. ")
    var_df = pd.read_table(args.SNP, sep = "\t", header = 0)
    gid = sorted(set(var_df["gene_id"]))
    ### split work
    if rank == 0:
        gid_num = len(gid)
        worker_tasks = {w:[] for w in range(size)}
        w_idx = 0
        for i in range(gid_num):
            worker_tasks[w_idx].append(i)
            w_idx = (w_idx + 1) % size
    else:
        worker_tasks = None
    worker_tasks = comm.bcast(worker_tasks, root = 0)
    if rank == 0:
        print("Split work. ")
    ###
    gene_file = open(args.output + "_tmp" + str(rank) + ".txt", "w")
    gene_file.write("gene_id\ttotal_read\talt1_read\talt2_read\tnon_read\tconflict\n")
    ### parallel
    for i in worker_tasks[rank]:
        vid, gene_read_count, alt1_read_count, alt2_read_count, non_read_count, cnft = variants_allele(args, gid, var_df, i, gene_dict, gstrand, bam_handle) # chr, gid, gstart, gend, var_num, total, alt1, alt2, non, conflict
        gene_file.write("{0:s}\t{1:d}\t{2:d}\t{3:d}\t{4:d}\t{5:d}\n".format(vid, var_num, gene_read_count, alt1_read_count, alt2_read_count, non_read_count, cnft))
    gene_file.close()
    df = pd.read_table(args.output + "_tmp" + str(rank) + ".txt", header = 0, sep = "\t")
    df = comm.gather(df, root = 0)
    if rank == 0:
        df = pd.concat(df, axis = 0)
        df = df.sort_values(by = ["contig", "gene_start", "gene_end"])
        df.to_csv(args.output + ".txt", sep = "\t", index = False)
        print("Write genotyping data successfully! ")
        check_rm(args, size)

##############
### Run it ###
##############

if __name__ == "__main__":
    main()
