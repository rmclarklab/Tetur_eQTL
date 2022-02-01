# Expression QTL (eQTL) of Tetranychus urticae
This is a repo for the data analysis pipeline of Tetranychus urticae eQTL.
eQTL analysis is an association analysis between genotype and gene expression.

# eQTL Project introduction
  To initiate eQTL project, we collected a total of 458 isogenic pools of the F3 T. urticae population. Briefly, a susceptible ROS-ITi and more resistant MR-VPi inbred strains are employed as the founder strains (F0). By crossing the two parental strains, we collected F1 female (diploid). And F1 female lay eggs without ferterlization developing into males (F2, haploid), which are back crossed to the susceptible ROS-ITi strain. For each F2 male backcross, all offsprings are collected to generate one single isogenic pool. 
  Because the recombination events happended during F1 reproducing F2 male, the F3 populations have different genotypic compositions which provided the foundamental basis for eQTL analysis. RNA was extracted from individual F3 isogenic pool, and RNA-seq was utilized for genotyping identification (See below for detail). 
  
# DNA-seq for variants calling
To call variants for the inbred ROS-ITi and MR-VPi strains, we mapped illumina DNA-seq against the reference genome (T. urticae, London strain). 
1. First, prepare index for the genome fasta file.
2. Then, map DNA-seq onto fasta genome using BWA.
3. Mark duplicated reads that are arised from PCR process.
4. Using the Best practice of GATK for variant calling. 
5. Filter SNP data for following processing.
6. Use the custom script to filter SNPs in homozygous genotype (given the parental strains are pure inbred).
7. (customized) To collect SNPs that are distinguishable between the two inbred parental strains (ROS-ITi vs. MR-VPi).

For each (no customized) step, you can refer to my other repo for detail. 

## In this step, we get ready the standard VCF file. SNPs in VCF file are used as diagnosable signal for the genotype call of RNA-seq of each F3 isogenic pool. 

# RNA-seq mapped against the reference genome. 
Data processing step by mapping RNA-seq onto the same reference genome (T. . urticae, London strain). 
1. First, prepare index file for the genome fasta file (no provided GTF, optional).
2. Map RNA-seq onto reference genome given the index folder. 

To pipeline the mapping process, see here. 

## In this step, we generated RNA-seq alignment file for each F3 isogenic pool in BAM format, which is required for the following analysis. 

# Call genotype composition based on RNA-seq alignment.
We developed a customized pipeline for the genotyping call of F3 isogenic pool in our study. 
Inputs:
- BAM file in coordinates sorted fashion and its index file.
- SNPs that are distinguishable between the two working inbred stains. 

Python with the following packages installed:
- pysam
- pandas
- mpi4py (to support multi-core running)

1. Processing SNPs genotype
2. Filter noises in SNP genotype call, and call genotype blocks for each isogenic pool
3. Combine all isogenic pool genotypic blocks, and generate recombination bin for following association analysis. 

## In this step, genotypes for individual F3 isogenic pools are generated. Genotype in recombination bin served as the representative genotype for following analysis. 

# Gene expression quantification.
For association analysis between genotype and gene expression, we need both data set be prepared. Aside from the file 



