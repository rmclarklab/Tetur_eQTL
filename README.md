# Expression QTL (eQTL) of <i>Tetranychus urticae</i> (a generalist spider-mite herbivore)
This is a repo for the data analysis pipeline of Tetranychus urticae eQTL project. </br>
eQTL is QTL explaining gene expression, can be identified via association analysis between genotype and gene expression.

## eQTL Project introduction
  To initiate eQTL project, we collected a total of 458 isogenic pools of the recombinant inbred lines (RIL). Briefly, a susceptible ROS-ITi (diploid mother, ♀) and more resistant MR-VPi (haploid father, ♂) inbred strains are employed as the founder strains (F0). By crossing the two parental strains, we collected F1 female (diploid). And F1 female lay eggs without ferterlization developing into males (F2, haploid), which are back crossed to the susceptible ROS-ITi strain. For each F2 male backcross, all offsprings are collected to generate one single isogenic pool. 
  Because the recombination events happended during F1 reproducing F2 male, the  RILs have different genotypic compositions which provided the foundamental basis for eQTL analysis. RNA was extracted from individual RIL isogenic populations, which are used for genotyping and phenotyping (phenotype data is gene expression level, See below for detail). 

## Procedure

- [DNA-seq for variants calling](#DNA-seq-for-variants-calling)
- [Map RNA-seq against the reference genome](#Map-RNA-seq-against-the-reference-genome)
- [Genotype call for RILs based on RNA-seq alignment](#Genotype-call-for-RILs-based-on-RNA-seq-alignment)
- [Update GFF3 file for the reference genome](#Update-GFF3-file-for-the-reference-genome)
- [Gene expression level quantification](#Gene-expression-level-quantification)
- [Association analysis between genotype and gene expression](Association-analysis-between-genotype-and-gene-expression)
- 

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
For Variants filtering, see [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants) also for hard-filtering. <br>

6. Pick SNPs that are distinguishable between the two inbred parental strains (ROS-ITi vs. MR-VPi). Output in tab-separated file
```bash
# Comparing filtered VCF files for ROS-IT and MR-VP, and pick genotype-calls different between them
vcf_compare.py -vcf1 ROS-IT.filtered.vcf.gz -vcf2 MR-VP.filtered.vcf.gz -R Tetranychus_urticae_2017.11.21.fasta -O variant_ROSIT.vs.MRVP
```
## Map RNA-seq against the reference genome
The three-chromosome reference genome was used, the same for DNA-seq mapping. <br>
We used two RNA-seq aligners, [STAR](https://github.com/alexdobin/STAR) and [GSNAP](https://github.com/juliangehring/GMAP-GSNAP), for RNA-seq mapping.

1. Generate indices for genome fasta file.
```bash
# STAR index generation
STAR --runMode genomeGenerate --runThreadN 30 --genomeDir STAR_index --genomeFastaFiles Tetranychus_urticae_2017.11.21.fasta --genomeSAindexNbases 12
# GSNAP index generation
gmap_build -D GSNAP_index -d Tetur Tetranychus_urticae_2017.11.21.fasta
```
2. Map RNA-seq onto the reference genome given the index folder. 
```bash
# STAR mapping, sort and index BAM alignment file
STAR --genomeDir STAR_index --runThreadN 20 --readFilesIn r1.fastq.gz r2.fastq.gz --twopassMode Basic --sjdbOverhang 99 --outFileNamePrefix sample_name. --readFilesCommand zcat --alignIntronMax 30000 --outSAMtype BAM Unsorted && samtools sort sample_name.Aligned.out.bam -o sample_name_sorted.bam -@ 8 && samtools index sample_name_sorted.bam 
```
  SNP-tolerant mapping using GSNAP require preparation of variant files
```bash
# make folder for SNP data
mkdir Tetur_SNP
# Format SNP allele information from vcf_compare.py output (see above) using SNP_prep.py
# SNP_prep.py [input] [output]
SNP_prep.py variant_ROSIT.vs.MRVP.txt SNP_allele
# Prepare SNP information for GSNAP
cat SNP_allele.txt | iit_store -o SNP_allele
mv SNP_allele.iit Tetur_SNP/Tetur_SNP.maps
# create a reference space index and compressed genome
snpindex -d Tetur_SNP -v SNP_allele -D . -V .
# prepare know splice sites 
cat gtf_file | gtf_splicesites > Tu.splicesites
cat gtf_file | gtf_introns > Tu.introns
cat Tu.splicesites | iit_store -o Tu_splicesites
cat Tu.introns | iit_store -o Tu_introns
# move all generated file into genome index db, then run GSNAP mapping
gsnap -d GSNAP_index -N 1 -D . --gunzip -s Tu_splicesites -v SNP_allele -t 20 -A sam r1.fastq.gz r2.fastq.gz | samtools sort -o sample_name.bam -O bam -@ 20 - && samtools index sample_name.bam
```
## Genotype call for RILs based on RNA-seq alignment
We developed a customized pipeline for genotyping purposes of RIL isogenic pools. 
Inputs:
  - BAM file in coordinates sorted fashion and with its index file;
  - SNPs information that are distinguishable between the two inbred stains (output of vcf_compare.py, see above).

1. Count allele-specific reads on SNP sites for each sample separately
```bash
# run genotype_allele.py to count allele-specifc reads on the SNP sites
# this is a multiple-core processing program, adjust core usage via "-n"
mpiexec -n 10 genotype_allele.py -V variant_ROSIT.vs.MRVP.txt -bam sample_name.bam -O sample_allele_count
```
After running for all samples, place all of them in the same folder (raw_count). <br>
2. Collect genotype information for all samples, and count the genotype frequency at each SNP site  
```bash
# run genotype_freq.py for genotype frequency, either in homozygous or heterozygous genotype
# this is a multiple-core processing program, adjust core usage via "-n"
mpiexec -n 10 genotype_freq.py -dir raw_count -O SNP_geno_freq
```
3. A backcrossing experimental design indicates a 1:1 ratio of heterozygous:homozygous genotype at each SNP site. We performed Chi-square goodness of fit test to filter bad SNP sites which doesn't fit the ratio (adjusted p < 0.01). <br>

    See chisq_bad.Rmd

4. For raw allele-specific read count, clean the dataset by dropping bad SNPs from last step
```bash
# run clean_count.R script to drop SNP rows in bad_SNP file
Rscript clean_count.R -raw sample_allele_count.txt -bad bad_SNPs.txt -O sample_allele_count.clean.txt
```
5. Call genotypic blocks based on allele-specific read count of good SNPs 
You need to set up the chromosomes of interested for genotype block assignment, and also provide chromosome length information in a tab-separated file. About how to provide chromosome length, see [here](https://biopython.org/docs/1.75/api/Bio.SeqIO.html). Example data set see under data folder.
```bash
# run genotype_block.py to call genotype blocks that arised from crossingover events 
genotype_block.py -chr chr.txt -chrLen chrlen.txt -C sample_allele_count.clean.txt -O sample_genotype_block
```
6. Visulization of genotype blocks on chromosome level
```bash
Rscript block_vis.R -geno sample_genotype_block.txt
```

## Update GFF3 file for the reference genome
To integrate all annotated gene information in the current three-chromosome reference genome, we added gene models from Orcae database and updated the current GFF3 file. <br>
To transfer gene models on fragmented scaffold genomes onto three-chromosome genome scale, see script ```restore_gff.py```. <br>
Combine the added gene models to the current GFF3 file and sort it using [gff3sort.pl](https://github.com/billzt/gff3sort). <br>
Transform from GFF3 to GTF format using script ```gff2gtf.py```.
```bash
# sort gff using gff3sort.pl
gff3sort.pl input.gff > output.gff
# compress gff
bgzip output.gff
# add index for compressed gff
tabix -p gff output.gff.gz
```

## Gene expression level quantification
By taking the updated GFF version, we run htseq-count on the alignment RNA-seq BAM files and output read count on gene basis.
```bash
# htseq-count command line
htseq-count -f bam -r pos -s reverse -t exon --nonunique none sample_name.bam Tetur.gtf > sample_name.txt
```

## Association analysis between genotype and gene expression
For each samples, its genotype blocks and gene expression data are available.
1. Merge htseq-count files for all samples into one single file
```bash
# make directory for htseq-count data of all samples
mkdir htseq_count_dir
# place all samples htseq-count output under the htseq_count_dir
mv *.txt htseq_count_dir
# Gene on row and sample on column
htseq_merge.py -countdir htseq_count_dir -O all_sample_count
```
2. Standarize read-count by performing library size normalization and quantile normalization. 

Library size normalization was performed using [DESeq2](https://github.com/mikelove/DESeq2). For running script, see ```htseq_norm.Rmd``` <br>

3. 

4. For any significant associations, recombination bins that are physically linked to each other are all passed the significance cutoff. To eliminate the issue arising from linkage disequilibruim (LD), we rebuild the linkage groups based on the bin genotype and then extracted the most significant association(s) between individual gene and its peak eQTL. 

First, we need to generate the a linkage measure for each bin to bin. 
Then, we developed a customized script to screening the output of MatrixeQTL. When one gene expression is associated with multiple recombination bins that belonged to one single linkage group, we only use the most significant association as the informative one. 



