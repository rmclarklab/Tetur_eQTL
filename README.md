# Expression QTL (eQTL) of <i>Tetranychus urticae</i> (a generalist spider-mite herbivor)
This is a repo for the data analysis pipeline of Tetranychus urticae eQTL project. </br>
eQTL is QTL explaining gene expression, can be identified via association analysis between genotype and gene expression.

## eQTL Project introduction
  To initiate eQTL project, we collected a total of 458 isogenic pools of the recombinant inbred lines (RIL). Briefly, a susceptible ROS-ITi (diploid mother, ♀) and more resistant MR-VPi (haploid father, ♂) inbred strains are employed as the founder strains (F0). By crossing the two parental strains, we collected F1 female (diploid). And F1 female lay eggs without ferterlization developing into males (F2, haploid), which are back crossed to the susceptible ROS-ITi strain. For each F2 male backcross, all offsprings are collected to generate one single isogenic pool. 
  Because the recombination events happended during F1 reproducing F2 male, the  RILs have different genotypic compositions which provided the foundamental basis for eQTL analysis. RNA was extracted from individual RIL isogenic populations, which are used for genotyping and phenotyping (phenotype data is gene expression level, See below for detail). 

## DNA-seq for variants calling
To call variants for the inbred ROS-ITi and MR-VPi strains, we mapped illumina DNA-seq against the three-chromosome reference genome (London strain, see [Wybouw, Kosterlitz, et al., 2019](https://academic.oup.com/genetics/article/211/4/1409/5931522)). <br>
GATK best practice for variants calling is refered [here](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows). <br>
1. First, prepare index for the genome fasta file;
```bash
# make directory for bwa index files
mkdir bwa_index
# change working directory to the folder
cd bwa_index
# generate index files using bwa index command
bwa index Tetranychus_urticae_2017.11.21.fasta
```
2. Then, map DNA-seq of ROS-ITi and MR-VPi onto the reference fasta genome using BWA;
```bash
# make directory for bwa mapping
mkdir bwa_map
# change working directory to the mapping folder
cd bwa_map
# run bwa mapping for ROS-ITi sample (paired-end DNA sequences)
bwa mem -t 20 -R "@RG\tID:20190412_8\tSM:ROS-ITi\tPL:Illumina\tLB:ROS-ITi" bwa_index/Tetranychus_urticae_2017.11.21.fasta r1.fastq.gz r2.fastq.gz | samtools view -Su - | samtools sort -@ 20 - -o ROS-ITi.BWA.bam
# run bwa mapping for MR-VPi sample (paired-end DNA sequences)
bwa mem -t 20 -R "@RG\tID:20190312\tSM:MR-VPi\tPL:Illumina\tLB:MR-VPi" bwa_index/Tetranychus_urticae_2017.11.21.fasta r1.fastq.gz r2.fastq.gz | samtools view -Su - | samtools sort -@ 20 - -o MR-VPi.BWA.bam
```
3. Mark duplicated reads that are arised from PCR using picard MarkDuplicate;
```bash
# mark duplicated reads in BWA mapping files
picard MarkDuplicates I=ROS-ITi.BWA.bam O=ROS-IT_duplicate.bam M=ROS-IT_metrics.txt && samtools index ROS-IT_duplicate.bam
picard MarkDuplicates I=MR-VPi.BWA.bam O=MR-VP_duplicate.bam M=MR-VP_metrics.txt && samtools index MR-VP_duplicate.bam
# left align insertion and deletion mappings (optional)
gatk LeftAlignIndels -R Tetranychus_urticae_2017.11.21.fasta -I ROS-IT_duplicate.bam -O ROS-IT_leftalign.bam
gatk LeftAlignIndels -R Tetranychus_urticae_2017.11.21.fasta -I MR-VP_duplicate.bam -O MR-VP_leftalign.bam
```
4. Run the Best practice of GATK pipeline for variant calling;
```bash
# run gatk HaplotypeCaller to generate g.vcf file
gatk HaplotypeCaller -R Tetranychus_urticae_2017.11.21.fasta -I ROS-IT_leftalign.bam -ERC GVCF -O ROS-IT.g.vcf.gz
gatk HaplotypeCaller -R Tetranychus_urticae_2017.11.21.fasta -I MR-VP_leftalign.bam -ERC GVCF -O MR-VP.g.vcf.gz
# run gatk GenotypeGVCFs to call variants in vcf file
gatk GenotypeGVCFs -R Tetranychus_urticae_2017.11.21.fasta -V ROS-IT.g.vcf.gz -O ROS-IT.vcf.gz
gatk GenotypeGVCFs -R Tetranychus_urticae_2017.11.21.fasta -V MR-VP.g.vcf.gz -O MR-VP.vcf.gz
```
5. Select variants data in unfiltered vcf file;
```bash
### run the following steps for ROS-ITi and MR-VPi, respectively
# select SNPs
gatk SelectVariants -R Tetranychus_urticae_2017.11.21.fasta -V input.vcf.gz -select-type-to-include SNP -O SNP.vcf.gz
# select INDELs (insertion and deletion, optional)
gatk SelectVariants -R Tetranychus_urticae_2017.11.21.fasta -V input.vcf.gz -select-type-to-include INDEL -O INDEL.vcf.gz
# sort vcf files
bcftools sort -o sorted.vcf.gz -O z SNP.vcf.gz
# add index for sorted vcf file
bcftools index -t sorted.vcf.gz
# Filter SNPs based on RMS mapping quality and genotype field information (run script vcf_pass.py)
vcf_pass.py -vcf sorted.vcf.gz -R Tetranychus_urticae_2017.11.21.fasta -O filtered.vcf.gz
```
For Variants filtering, [see](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants) also for hard-filtering. 

7. Collect SNPs that are distinguishable between the two inbred parental strains (ROS-ITi vs. MR-VPi).
```bash

```

## Update GFF3 file for the reference genome
To integrate all annotated genes in the current reference genome, we provided a newer version of GFF3 annotation file for the working reference genome. 


In this step, we prepared the standard VCF file. SNPs in VCF file are used as diagnosable signal for the genotype call of RNA-seq of each RIL isogenic pools. 

## RNA-seq mapped against the reference genome. 
Data processing step by mapping RNA-seq onto the three-chromosome reference genome (same as we used for DNA-seq alignment). 
1. First, prepare index file for the genome fasta file using STAR.
2. Map RNA-seq onto reference genome given the index folder using STAR mapping. 

To pipeline the mapping process, see here. 

#### In this step, we generated RNA-seq alignment file (BAM) for each RIL isogenic pool, which is required for the following analysis. 

## Call genotype composition based on RNA-seq alignment.
We developed a customized pipeline for the genotyping call of RIL isogenic pool in our study. 
[[Inputs]]
- BAM file in coordinates sorted fashion and its index file.
- SNPs that are distinguishable between the two working inbred stains. 

Python with the following packages installed:
- pysam
- pandas
- mpi4py (to support multi-core running)

1. Processing SNPs genotype
2. Filter noises in SNP genotype call, and call genotype blocks for each isogenic pool
3. Combine all isogenic pool genotypic blocks, and generate recombination bin for following association analysis. 

#### In this step, genotypes for individual F3 isogenic pools are generated. Genotype in recombination bin served as the representative genotype for following analysis. 

## Gene expression quantification.
Aside from the genotype data, we need to generate gene expression data for association analysis between them. 
Here, we still use the RNA-seq alignment file in BAM format for quantify gene expression. 
Using htseq-count to count expression on gene-basis. 

## Association analysis between genotype and gene expression.
We used MatrixeQTL for the association analysis between genotype and gene expression. 
Inputs:
- Genotype on recombination bins;
- Gene expression of all isogenic pools.
Command line:

## Significant association extraction
For any significant associations, recombination bins that are physically linked to each other are all passed the significance cutoff. To eliminate the issue arising from linkage disequilibruim (LD), we rebuild the linkage groups based on the bin genotype and then extracted the most significant association(s) between individual gene and its peak eQTL. 

1. First, we need to generate the a linkage measure for each bin to bin. 
2. Then, we developed a customized script to screening the output of MatrixeQTL. When one gene expression is associated with multiple recombination bins that belonged to one single linkage group, we only use the most significant association as the informative one. 

#### In this step, we parsed the output of MatrixeQTL and only collected the most significant association for following check. 




