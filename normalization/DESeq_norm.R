###########################################################
### This script is for differential expression analysis
### with DESeq2. Accept read count data file and 
### column information. 
###########################################################

suppressMessages(library("DESeq2"))
suppressMessages(library("argparse"))
suppressMessages(library("stringr"))

parse_arg <- function() {
    # accept arguments from terminal
    parser <- ArgumentParser(description = "Given sample read count files, return normalized read count data in one combined file.")
    parser$add_argument("-count", "--count", nargs = "+", help = "read count files. ")
    parser$add_argument("-O", "--output", help = "output prefix. ")
    parser$parse_args()
}

###
argv <- parse_arg()
count_files <- argv$count
out <- argv$output

### combine all read count files in one file (column name using file name)
df_list <- list()
for (f in count_files) {
    file_base <- basename(f)
    pre <- strsplit(file_base, ".txt")[[1]]
    cnt_df <- read.table(f, header = F, sep = "\t", row.names = 1)
    colnames(cnt_df) <- pre
    df_list[[pre]] <- cnt_df
}

df_all <- do.call(cbind, df_list)

### generate sample file based on read count data 
samples <- colnames(df_all)
sample_df <- data.frame(sample = samples)
sample_df$genotype <- str_sub(sample_df$sample, start = 1, end = -3)
# sample_df$replicate <- str_sub(sample_df$sample, start = -1)

### use DESeq2 to normalizing read count 
dds <- DESeqDataSetFromMatrix(countData = df_all, colData = sample_df, design = ~ 1)
dds <- estimateSizeFactors(dds)
norm_cts <- counts(dds, normalized = T)

# to transform the count data to the log2 scale which minimizes differences between samples for rows with sample counts, and normalizes with respect to library size (rlog). 
# rld <- rlog(dds)
# norm_cts <- as.data.frame(assay(rld))

write.table(norm_cts, file = paste0(out, ".txt"), sep = "\t", quote = F, row.names = T)
