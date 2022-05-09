
import argparse
import os
import pandas as pd

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description=" \n\
        Usage: python allele_cis.py -countdir <cnt> -genodir <geno> -gene <gene> -O <prefix>")
    parser.add_argument("-countdir", "--countdir", help = "directory for count files. ")
    parser.add_argument("-genodir", "--genodir", help="directory for genotype files. ")
    parser.add_argument("-gene", "--gene", help = "the gene information (including start and end position) file. ")
    parser.add_argument("-O", "--output", help="output table with allele-specific read count on gene-basis. ")
    args = parser.parse_args()
    return args

def check_fs(args):
    '''check the files under countdir and genodir are identical'''
    countdir = args.countdir
    genodir = args.genodir
    countf = sorted([f for f in os.listdir(countdir) if f.endswith(".txt")])
    genof = sorted([f for f in os.listdir(genodir) if f.endswith(".txt")])
    if countf == genof:
        print("The input count files and genotype files are matched!")
        return True
    else:
        print("The input count files and genotype files are not matched! Please check! ")
        return False

def search_heter(args):
    ###
    genodir = args.genodir
    countdir = args.countdir
    gene = args.gene
    output = args.output
    ###
    genof = sorted([f for f in os.listdir(genodir) if f.endswith(".txt")])
    gloc = pd.read_table(gene, sep = "\t", header = 0)
    gloc["snp"] = gloc["contig"] + ":" + gloc["start"].astype(str)
    ###
    sample_het = list()
    for g in genof:
        pre = g.split(".txt")[0]
        gdf = pd.read_table(os.path.join(genodir, g), sep = "\t", header = 0)
        heteloc = gdf[gdf["genotype"] == "hete"]
        chrhet = {row.chromosome:(row.start,row.end) for row in heteloc.itertuples()}
        hetgid_list = list()
        for c in chrhet:
            sta = chrhet[c][0]
            end = chrhet[c][-1]
            hetgid = list(gloc[(gloc["contig"] == c) & (gloc["start"] >= sta) & (gloc["end"] <= end)]["gene_id"])
            hetgid_list = hetgid_list + hetgid
        ###
        cdf = pd.read_table(os.path.join(countdir, g), sep = "\t", header = 0)
        hetcnt = cdf[cdf["gene_id"].isin(hetgid_list)]
        hetcnt = hetcnt[((hetcnt["non_read"] + hetcnt["conflict"])/hetcnt["total_read"] < 0.05) & (hetcnt["total_read"] > 10)]
        hetcnt = hetcnt[["gene_id", "alt1_read", "alt2_read"]]
        hetcnt["total"] = hetcnt["alt1_read"] + hetcnt["alt2_read"]
        hetcnt = hetcnt[["gene_id", "alt1_read", "total"]]
        hetcnt["id"] = pre
        hetdf = hetcnt.merge(gloc[["gene_id", "snp"]], on = "gene_id", how = "left")
        sample_het.append(hetdf)
    het_gene = pd.concat(sample_het, axis = 0)
    het_gene = het_gene[["gene_id", "id", "snp", "alt1_read", "total"]]
    het_gene = het_gene.rename(columns = {"gene_id":"gene", "alt1_read":"ref"})
    het_gene.to_csv(output + ".txt", sep = "\t", index = False)

def main():
    args = args_parser()
    if check_fs(args):
        search_heter(args)

##############
### Run it ###
##############

if __name__ == "__main__":
    main()