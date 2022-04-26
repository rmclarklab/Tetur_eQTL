import os 
import pandas as pd 
import pysam 
from Bio import SeqIO

sanger_fa = "Orcae_01252019/Tetranychus_urticae.main_genome_200909.scaffolds.fasta"
pseudo_chr_fa = "Clark_2019/Tetranychus_urticae_2017.11.21.fasta"

sanger_seq = SeqIO.parse(sanger_fa, "fasta")
pseudo_chr_seq = SeqIO.parse(pseudo_chr_fa, "fasta")

sanger_dict = SeqIO.to_dict(sanger_seq)
pseudo_chr_dict = SeqIO.to_dict(pseudo_chr_seq)

### break out the scaffold fasta file given breakpoint information
sangerbreak = pd.read_table("data/sangerbreaks.txt", sep = "\t", header = 0)
scaf_list = sorted(set(sangerbreak["scaffold_old"]))

scaf_new_dict_sub = dict()
for s in scaf_list:
    sub_scaf = sangerbreak[sangerbreak["scaffold_old"] == s]
    seq_str = str(sanger_dict[s].seq)
    seq_len = len(seq_str)
    breakpoint_list = list()
    for row in sub_scaf.itertuples():
        rstart = row.start
        rend = row.end
        if rend != "-":
            breakpoint_list.append(int(rend))
        else:
            breakpoint_list.append(seq_len)
    b0 = 0
    i = 1
    for b in breakpoint_list:
        sub_seq = seq_str[b0:b]
        b0 = b
        scaf_new_dict_sub[s+"."+str(i)] = sub_seq
        i += 1

scaf_new_fw = open("data/Tetranychus_urticae.main_genome_200909.scaffolds.break.fasta", "w")

for entry in sanger_dict:
    sid = sanger_dict[entry].id
    if sid not in scaf_list:
        sseq = str(sanger_dict[entry].seq)
        scaf_new_fw.write(f">{sid}\n{sseq}\n")
    else:
        bb_key = [bb for bb in scaf_new_dict_sub if sid in bb]
        for bb in bb_key:
            sseq = scaf_new_dict_sub[bb]
            scaf_new_fw.write(f">{bb}\n{sseq}\n")

scaf_new_fw.close()

### put scaffold back onto chromosome to test the number all correct
def reverse(seq):
    """Returns a reversed string"""
    return seq[::-1]

def complement(seq):
    """Returns a complement DNA sequence"""
    complement_dict = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'N':'N'}
    seq_list = list(seq)
    seq_list = [complement_dict[base] for base in seq_list]
    return ''.join(seq_list)

scaf_new = "data/Tetranychus_urticae.main_genome_200909.scaffolds.break.fasta"
scaf_new_seq = SeqIO.parse(scaf_new, "fasta")
scaf_new_dict = SeqIO.to_dict(scaf_new_seq)

scaf_order = pd.read_table("data/scaffold_order.txt", sep = "\t", header = 0)
chr_set = sorted(set(scaf_order["chromosome"]))
scaf_chr_set = sorted(set(scaf_order["scaffold"]))
chr_new_fw = open("data/chromosome_linked.fasta", "w")

for c in chr_set:
    scaf_c = scaf_order[scaf_order["chromosome"] == c]
    chr_scaffold_link = list()
    for row in scaf_c.itertuples():
        rscaffold = row.scaffold
        rdirection = row.direction
        scaf_c_seq = str(scaf_new_dict[rscaffold].seq)
        if rdirection == "-":
            scaf_c_seq = complement(reverse(scaf_c_seq))
        chr_scaffold_link.append(scaf_c_seq)
    chr_linked = "".join(chr_scaffold_link)
    chr_new_fw.write(f">chromosome_{c}\n{chr_linked}\n")

for entry in scaf_new_dict:
    eid = scaf_new_dict[entry].id
    eseq = scaf_new_dict[entry].seq
    if eid not in scaf_chr_set:
        chr_new_fw.write(f">{eid}\n{eseq}\n")

chr_new_fw.close()

# compare to andre's chromosome 

mychr = "data/chromosome_linked.fasta"
mychr_seq = SeqIO.parse(mychr, "fasta")
mychr_dict = SeqIO.to_dict(mychr_seq)

for entry in mychr_dict:
    mseq = str(mychr_dict[entry].seq)
    aseq = str(pseudo_chr_dict[entry].seq)
    if mseq == aseq:
        print(1)
    else:
        print(0)
# exactly the same

### add scaffold to chromosome mapping relationship

scaf_new = "data/Tetranychus_urticae.main_genome_200909.scaffolds.break.fasta"
scaf_new_seq = SeqIO.parse(scaf_new, "fasta")
scaf_new_dict = SeqIO.to_dict(scaf_new_seq)

scaf_order = pd.read_table("data/scaffold_order.txt", sep = "\t", header = 0)
chr_set = sorted(set(scaf_order["chromosome"]))
scaf_chr_set = sorted(set(scaf_order["scaffold"]))

scaffold_range_dict = dict()

for c in chr_set:
    chri = 0
    scaf_c = scaf_order[scaf_order["chromosome"] == c]
    for row in scaf_c.itertuples():
        rscaffold = row.scaffold
        rscaffold_seq = str(scaf_new_dict[rscaffold].seq)
        rscaffold_len = len(rscaffold_seq)
        scaffold_range_dict[rscaffold] = rscaffold_len
        scaffold_range_dict[rscaffold + ".start"] = chri + 1
        chri = chri + rscaffold_len
        scaffold_range_dict[rscaffold + ".end"] = chri

scaffold_len_list = [scaffold_range_dict[row.scaffold] for row in scaf_order.itertuples()]
scaffold_start = [scaffold_range_dict[row.scaffold + ".start"] for row in scaf_order.itertuples()]
scaffold_end = [scaffold_range_dict[row.scaffold + ".end"] for row in scaf_order.itertuples()]

scaf_order["scaffold_len"] = scaffold_len_list
scaf_order["start"] = scaffold_start
scaf_order["end"] = scaffold_end

scaf_order.to_csv("data/scaffold_order_new.txt", sep = "\t", index = False)

### read Orcae GFF3 and Clark GFF3 for unique geneid in each dataset
gff3_orcae = pd.read_table("Orcae_01252019/scaffold.gff3", comment = "#", header = None)
gff3_clark = pd.read_table("Clark_2019/Tetranychus_urticae_2018.06.01.gff3", comment = "#", header = None)

def gff_geneid(gff):
    gff = gff.rename(columns = {0:"tig", 1:"source", 2:"feature", 3:"start", 4:"end", 5:"score", 6:"strand", 7:"phase", 8:"attributes"})
    gene_gff = gff[gff["feature"] == "gene"]
    geneid_list = list(gene_gff.attributes.str.split(";", expand = True)[[0]][0])
    return geneid_list

orcae_geneid = gff_geneid(gff3_orcae)
clark_geneid = gff_geneid(gff3_clark)

orcae_geneid = [o.split("=")[-1] for o in orcae_geneid]
clark_geneid = [c.split(":")[-1] for c in clark_geneid]

orcae_unique = sorted(set(orcae_geneid).difference(set(clark_geneid))) # 45
clark_unique = sorted(set(clark_geneid).difference(set(orcae_geneid))) # 18

list_df = list()
for o in orcae_unique:
    gff3_sub = gff3_orcae[gff3_orcae[8].str.contains(o)]
    list_df.append(gff3_sub)

gff3_orcae_use = pd.concat(list_df, axis = 0)

gff3_orcae_use.to_csv("data/Orcae_unique.gff3", sep = "\t", index = False)

### shift the location from scaffold to chromosome
gff3_orcae_use = pd.read_table("data/Orcae_unique.gff3", sep = "\t", header = 0)
gff3_orcae_use = gff3_orcae_use.rename(columns = {"0":"tig", "1":"source", "2":"feature", "3":"start", "4":"end", "5":"score", "6":"strand", "7":"phase", "8":"attributes"})
coord = pd.read_table("data/scaffold2chr_coordinate.txt", sep = "\t", header = 0)

range_dict = dict()
for row in coord.itertuples():
    scaff_post = row.scaffold_break
    scaff_direction = row.direction
    if scaff_direction == "+":
        range_dict[scaff_post] = [range(row.scaffold_start, row.scaffold_end, 1), range(row.start, row.end, 1)]
    elif scaff_direction == "-":
        range_dict[scaff_post] = [range(row.scaffold_end, row.scaffold_start, -1), range(row.start, row.end, 1)]
    else:
        print("strand not indicated. ")

row_list = list()
strand_dict = {"+":"-", "-":"+", ".":"."}
for index, row in gff3_orcae_use.iterrows():
    rtig = row.tig
    rstart = row.start
    rend = row.end
    rstrand = row.strand
    scaf_search_sub = coord[(coord["scaffold_before_break"] == rtig) & (coord["scaffold_start"] < rstart) & (coord["scaffold_end"] > rend)]
    if len(scaf_search_sub) > 0:
        scaf_search = scaf_search_sub["scaffold_break"].item()
        chr_search = scaf_search_sub["chromosome"].item()
        sstart = range_dict[scaf_search][0].index(rstart)
        cstart = list(range_dict[scaf_search][1])[sstart]
        send = range_dict[scaf_search][0].index(rend)
        cend = list(range_dict[scaf_search][1])[send]
        row.tig = "chromosome_" + str(chr_search)
        row.start = cstart
        row.end = cend
        if row.start > row.end:
            row.strand = strand_dict[rstrand]
            row.start = cend
            row.end = cstart
    row_dict = row.to_dict()
    row_list.append(row_dict)

df_new = pd.DataFrame(row_list)
df_new.to_csv("data/Orcae_unique_pos_chr.gff3", sep = "\t", index = False)

### format gff3 extra gene information

gff_new = pd.read_table("data/Orcae_unique_pos_chr_manual.gff3", sep = "\t", header = 0)
parent_dict = {"CDS":"mRNA", "exon":"mRNA", "mRNA":"gene", "three_prime_UTR":"mRNA", "five_prime_UTR":"mRNA"}

gene_list_record = list()
gff3_format = open("data/Orcae_unique_pos_chr_format.gff3", "w")

for index, row in gff_new.iterrows():
    rfeature = row.feature
    rattributes = row.attributes
    if rfeature != "gene":
        id = rattributes.split("ID=")[-1].split(";")[0]
        parent = rattributes.split("Parent=")[-1].split(";")[0]
        id_label = "ID=" + rfeature +  ":" + id + ";"
        name_label = "Parent=" + parent_dict[rfeature] + ":" + parent + ";"
        parent_left = ";".join(rattributes.split("Parent=")[-1].split(";")[1:])
        new_attributes = id_label + name_label + parent_left
    else:
        id = rattributes.split("ID=")[-1].split(";")[0]
        id_left = ";".join(rattributes.split("ID=")[-1].split(";")[1:])
        id_label = "ID=" + rfeature + ":" + id + ";"
        new_attributes = id_label + id_left
        if id not in gene_list_record:
            gene_list_record.append(id)
            gff3_format.write("###\n")
    row.attributes = new_attributes
    row_l = [str(r) for r in row.tolist()]
    row_w = "\t".join(row_l)
    gff3_format.write(f"{row_w}\n")

gff3_format.close()

